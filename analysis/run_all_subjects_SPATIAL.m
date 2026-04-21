% run_all_subjects_NO_GLOBAL.m
% Same as run_all_subjects_GLOBAL but WITHOUT global luminance term
% This is a model comparison: local patches only vs. local + global

clear all
close all

% ====================================================================
% UPDATE THESE PATHS:
% ====================================================================
code_dir = '/Users/jaspreetdodd/Desktop/local analysis/code';
data_dir = '/Users/jaspreetdodd/Desktop/local analysis/data';
results_dir = fullfile(data_dir, 'results_no_global');  % Different folder!
% ====================================================================

addpath(code_dir);
cd(data_dir);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

subjects = 1;

fprintf('\n=========================================\n');
fprintf('ANALYSIS WITHOUT GLOBAL LUMINANCE TERM\n');
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
    
    tic;
    
    clearvars -except subj_num subjects code_dir data_dir results_dir summary_data
    
    %% Load data
    data_file = sprintf('subj%02d_retinalcoor.mat', subj_num);
    if ~exist(data_file, 'file')
        fprintf('WARNING: %s not found! Skipping...\n', data_file);
        continue;
    end
    load(data_file)
    
    op = 2; % Still use pupil model

    %% Set up data - WITHOUT GLOBAL LUMINANCE
    d = 0; 
    ig = 120;

    X = []; Y = [];
    for i = 1:length(matrix),
        % KEY CHANGE: Remove global_screen term, only keep intercept + patches
        X = [X; [ones(size(global_screen{i}(ig+1:end-d))) matrix{i}(ig+1:end-d,:)]];
        
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
    for i = 2:size(X,2),  % Start at 2 because column 1 is intercept
        nanx = isnan(X(:,i));
        X(nanx,i) = interp1(t(~nanx), X(~nanx,i), t(nanx), 'linear', 'extrap');
        [fs, fp] = pupil_create_command(X(:,i)', arrousal', SCinput', taup, taus, dt);
        P(:,i) = pupil_model(fs, fp, dt);
        
        if mod(i, 500) == 0
            fprintf('  Converting luminance: %1.1f%%\n', 100*i/size(X,2));
        end
    end

    %% Cross-validate regularized linear regression model
    % THIS IS THE "regularization hyperparameter optimization" Gunnar mentioned
    Lambda = logspace(-5,5,30);
    fprintf('  Running cross-validation over %d lambda values...\n', length(Lambda));
    
    CVMdl = fitrlinear(P,Y,'ObservationsIn','rows','KFold',5,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    mse = kfoldLoss(CVMdl);

    % Find optimal lambda (the optimization step!)
    ind = find(mse == min(mse));
    fprintf('  Optimal lambda selected: %.6f (index %d/%d)\n', Lambda(ind), ind, length(Lambda));
    
    % Train final model with optimal lambda
    Mdl = fitrlinear(P,Y,'ObservationsIn','rows','Lambda',Lambda(ind),...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');

    % Compute R^2
    Yp = predict(Mdl,P);
    SSR = nansum((Yp - Y).^2);
    TSS = nansum(((Y - nanmean(Y)).^2));
    Rsquared = 1 - SSR/TSS;
    
    fprintf('  R² = %.4f (without global), λ = %.6f\n', Rsquared, Lambda(ind));

    %% Save results
    % Extract spatial weights (now starts at index 2, not 3!)
    spatial_weights = Mdl.Beta(2:end);  % Skip only intercept (no global term)
    
    % Reshape to spatial map
    weight_map = reshape(spatial_weights, metadata.nY, metadata.nX);
    
    % Compute summary statistics
    n_sig_patches = sum(abs(spatial_weights) > 0.001);
    
    % Store in summary arrays
    summary_data.subject(end+1) = subj_num;
    summary_data.r2(end+1) = Rsquared;
    summary_data.optimal_lambda(end+1) = Lambda(ind);
    summary_data.min_mse(end+1) = min(mse);
    summary_data.n_significant_patches(end+1) = n_sig_patches;
    
    %% Save compact results
    save(fullfile(results_dir, sprintf('subj%02d_results_no_global.mat', subj_num)), ...
        'weight_map', 'Rsquared', 'Lambda', 'ind', 'mse', 'metadata', 'Mdl', '-v7.3');
    
    %% Save weight map figure
    fig = figure('Position', [100 100 800 600], 'Visible', 'off');
    imagesc(weight_map)
    colormap(parula)
    c = colorbar;
    caxis(max(abs(c.Limits))*[-1 1])
    title(sprintf('Subject %d (NO GLOBAL): R² = %.4f', subj_num, Rsquared))
    xlabel('Horizontal Patch Index')
    ylabel('Vertical Patch Index')
    
    saveas(fig, fullfile(results_dir, sprintf('subj%02d_weights_no_global.png', subj_num)));
    close(fig);
    
    elapsed = toc;
    fprintf('  ✓ Completed in %.1f minutes\n', elapsed/60);
    
end

%% Save summary
fprintf('\n========================================\n');
fprintf('CREATING SUMMARY STATISTICS\n');
fprintf('========================================\n');

summary_table = table(...
    summary_data.subject(:), ...
    summary_data.r2(:), ...
    summary_data.optimal_lambda(:), ...
    summary_data.min_mse(:), ...
    summary_data.n_significant_patches(:), ...
    'VariableNames', {'Subject', 'R2', 'OptimalLambda', 'MinMSE', 'NumSignificantPatches'});

writetable(summary_table, fullfile(results_dir, 'summary_NO_GLOBAL.csv'));
save(fullfile(results_dir, 'summary_NO_GLOBAL.mat'), 'summary_table', 'summary_data');

fprintf('\n✓ Summary saved to: %s\n', fullfile(results_dir, 'summary_NO_GLOBAL.csv'));

%% Display statistics
fprintf('\n========================================\n');
fprintf('FINAL SUMMARY - NO GLOBAL TERM (N = %d)\n', length(summary_data.r2));
fprintf('========================================\n');
fprintf('R² Statistics:\n');
fprintf('  Mean ± SD: %.4f ± %.4f\n', mean(summary_data.r2), std(summary_data.r2));
fprintf('  Range: [%.4f, %.4f]\n', min(summary_data.r2), max(summary_data.r2));
fprintf('  Median: %.4f\n', median(summary_data.r2));
fprintf('\n========================================\n');
fprintf('Compare this to the WITH GLOBAL results!\n');
fprintf('Expected: R² should be slightly lower without global term\n');
fprintf('========================================\n');
