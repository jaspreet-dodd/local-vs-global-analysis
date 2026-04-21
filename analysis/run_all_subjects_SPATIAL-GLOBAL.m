% run_all_subjects_GLOBAL.m
% Runs analysis for all 51 subjects

clear all
close all

% ====================================================================
% UPDATE THESE PATHS:
% ====================================================================
code_dir = '/Users/jaspreetdodd/Desktop/local analysis/code';
data_dir = '/Users/jaspreetdodd/Desktop/local analysis/data';
results_dir = fullfile(data_dir, 'results');
% ====================================================================

addpath(code_dir);
cd(data_dir);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

subjects = 1:51;  % ALL subjects!

fprintf('\n=========================================\n');
fprintf('OPTIMIZED ANALYSIS: Saving only essential data\n');
fprintf('Processing %d subjects\n', length(subjects));
fprintf('=========================================\n\n');

% Pre-allocate summary arrays
summary_data = struct();
summary_data.subject = [];
summary_data.r2 = [];
summary_data.optimal_lambda = [];
summary_data.min_mse = [];
summary_data.n_significant_patches = [];

for subj_num = subjects
    fprintf('\n========================================\n');
    fprintf('Processing Subject %d of %d\n', subj_num, length(subjects));
    fprintf('========================================\n');
    
    tic; % Time the processing
    
    % Clear workspace but keep our loop variables
    clearvars -except subj_num subjects code_dir data_dir results_dir summary_data
    
    %% Load data
    data_file = sprintf('subj%02d_retinalcoor.mat', subj_num);
    if ~exist(data_file, 'file')
        fprintf('WARNING: %s not found! Skipping...\n', data_file);
        continue;
    end
    load(data_file)
    
    op = 2; % pupil-simulated model

    %% Set up data
    d = 0; 
    ig = 120;

    X = []; Y = [];
    for i = 1:length(matrix),
        X = [X; [ones(size(global_screen{i}(ig+1:end-d))) global_screen{i}(ig+1:end-d) matrix{i}(ig+1:end-d,:)]];
        
        y = pupil_mm{i};
        nanx = isnan(y);
        t = 1:numel(y);
        y(nanx) = interp1(t(~nanx), y(~nanx), t(nanx), 'linear', 'extrap');
        Yd = detrend(y);
        ff = normpdf(-10:10,0,.1); ff = ff./sum(ff);
        Yd = conv(Yd, ff,'same');
        Y = [Y; Yd((ig+1+d):end)];
    end

    %% Convert luminance to pupil size
    dt = 0.0167;
    taus = 0.64;
    taup = 0.15;

    P = ones(size(X));
    arrousal = zeros(size(X,1),1);
    SCinput = zeros(size(X,1),1);
    
    t = 1:numel(X(:,1));
    for i = 2:size(X,2),
        nanx = isnan(X(:,i));
        X(nanx,i) = interp1(t(~nanx), X(~nanx,i), t(nanx), 'linear', 'extrap');
        [fs, fp] = pupil_create_command(X(:,i)', arrousal', SCinput', taup, taus, dt);
        P(:,i) = pupil_model(fs, fp, dt);
        
        if mod(i, 500) == 0
            fprintf('  Converting luminance: %1.1f%%\n', 100*i/size(X,2));
        end
    end

    %% Cross-validate regularized linear regression model
    Lambda = logspace(-5,5,30);
    fprintf('  Running cross-validation...\n');
    
    CVMdl = fitrlinear(P,Y,'ObservationsIn','rows','KFold',5,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    mse = kfoldLoss(CVMdl);

    ind = find(mse == min(mse));
    
    Mdl = fitrlinear(P,Y,'ObservationsIn','rows','Lambda',Lambda(ind),...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');

    % Compute R^2
    Yp = predict(Mdl,P);
    SSR = nansum((Yp - Y).^2);
    TSS = nansum(((Y - nanmean(Y)).^2));
    Rsquared = 1 - SSR/TSS;
    
    fprintf('  R² = %.4f, λ = %.6f\n', Rsquared, Lambda(ind));

    %% ================================================================
    %  SAVE ONLY ESSENTIAL DATA (NOT THE HUGE MATRICES)
    %  ================================================================
    
    % Extract spatial weights (the main result!)
    spatial_weights = Mdl.Beta(3:end);  % Skip intercept and global
    
    % Reshape to spatial map
    weight_map = reshape(spatial_weights, metadata.nY, metadata.nX);
    
    % Compute summary statistics
    n_sig_patches = sum(abs(spatial_weights) > 0.001);  % Count non-zero weights
    
    % Store in summary arrays
    summary_data.subject(end+1) = subj_num;
    summary_data.r2(end+1) = Rsquared;
    summary_data.optimal_lambda(end+1) = Lambda(ind);
    summary_data.min_mse(end+1) = min(mse);
    summary_data.n_significant_patches(end+1) = n_sig_patches;
    
    %% Save COMPACT results file
    % Only save what you actually need for figures!
    save(fullfile(results_dir, sprintf('subj%02d_results_compact.mat', subj_num)), ...
        'weight_map', ...           % The spatial weights (49x62 matrix)
        'Rsquared', ...             % R² value
        'Lambda', ...               % All lambda values tested
        'ind', ...                  % Index of optimal lambda
        'mse', ...                  % Cross-validation MSE curve
        'metadata', ...             % Grid dimensions etc.
        'Mdl', ...                  % Model object (for predictions if needed)
        '-v7.3');
    
    % Note: We're NOT saving X, Y, P, Yp (the huge data matrices)
    % This reduces file size from ~1GB to ~5-10MB per subject!
    
    %% Create and save weight map figure
    fig = figure('Position', [100 100 800 600], 'Visible', 'off');
    imagesc(weight_map)
    colormap(parula)
    c = colorbar;
    caxis(max(abs(c.Limits))*[-1 1])
    title(sprintf('Subject %d: R² = %.4f', subj_num, Rsquared))
    xlabel('Horizontal Patch Index')
    ylabel('Vertical Patch Index')
    
    saveas(fig, fullfile(results_dir, sprintf('subj%02d_weights.png', subj_num)));
    close(fig);
    
    elapsed = toc;
    fprintf('  ✓ Completed in %.1f minutes\n', elapsed/60);
    
end

%% Save overall summary
fprintf('\n========================================\n');
fprintf('CREATING SUMMARY STATISTICS\n');
fprintf('========================================\n');

% Convert to table
summary_table = table(...
    summary_data.subject(:), ...
    summary_data.r2(:), ...
    summary_data.optimal_lambda(:), ...
    summary_data.min_mse(:), ...
    summary_data.n_significant_patches(:), ...
    'VariableNames', {'Subject', 'R2', 'OptimalLambda', 'MinMSE', 'NumSignificantPatches'});

% Save as CSV
writetable(summary_table, fullfile(results_dir, 'ALL_SUBJECTS_summary.csv'));

% Also save as MAT for easy loading
save(fullfile(results_dir, 'ALL_SUBJECTS_summary.mat'), 'summary_table', 'summary_data');

fprintf('\n✓ Summary saved to: %s\n', fullfile(results_dir, 'ALL_SUBJECTS_summary.csv'));

%% Display summary statistics
fprintf('\n========================================\n');
fprintf('FINAL SUMMARY (N = %d subjects)\n', length(summary_data.r2));
fprintf('========================================\n');
fprintf('R² Statistics:\n');
fprintf('  Mean ± SD: %.4f ± %.4f\n', mean(summary_data.r2), std(summary_data.r2));
fprintf('  Range: [%.4f, %.4f]\n', min(summary_data.r2), max(summary_data.r2));
fprintf('  Median: %.4f\n', median(summary_data.r2));
fprintf('\nLambda Statistics:\n');
fprintf('  Mean: %.6f\n', mean(summary_data.optimal_lambda));
fprintf('  Median: %.6f\n', median(summary_data.optimal_lambda));
fprintf('\n========================================\n');
fprintf('ALL PROCESSING COMPLETE!\n');
fprintf('Results saved in: %s\n', results_dir);
fprintf('\nTotal storage used: ~%.1f MB (instead of ~%.1f GB)\n', ...
    length(summary_data.r2) * 5, length(summary_data.r2));
fprintf('========================================\n');
