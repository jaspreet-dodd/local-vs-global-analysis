% run_global_only_model.m
% Fits a simple global-luminance-only model for each subject
% (intercept + mean screen luminance, no spatial patches)
% Used to compare against spatial-only model R²
%
% Run this AFTER main analysis is complete.
% Results saved to results_dir/results_global_only/

clear; close all;

% ====================================================================
% UPDATE THESE PATHS
% ====================================================================
code_dir         = '/Users/jaspreetdodd/Desktop/local analysis/code/actually used';
data_dir         = '/Users/jaspreetdodd/Desktop/local analysis/data';
results_dir_nog  = '/Users/jaspreetdodd/Desktop/results/results_no_global';
results_dir_out  = '/Users/jaspreetdodd/Desktop/results/results_global_only';
% ====================================================================

addpath(code_dir);
if ~exist(results_dir_out, 'dir'), mkdir(results_dir_out); end

good_subjects = [1:17 19:51];
subjects      = good_subjects;

fprintf('=============================================\n');
fprintf('GLOBAL-ONLY MODEL\n');
fprintf('Intercept + mean screen luminance only\n');
fprintf('N = %d subjects\n', length(subjects));
fprintf('=============================================\n\n');

r2_global_only = nan(1, length(subjects));

for idx = 1:length(subjects)
    s  = subjects(idx);

    % Load retinalcoor data (need global_screen and pupil_mm)
    f_data = fullfile(data_dir, sprintf('subj%02d_retinalcoor.mat', s));
    if ~exist(f_data, 'file')
        fprintf('  WARNING: retinalcoor not found for subject %d, skipping\n', s);
        continue;
    end
    load(f_data, 'global_screen', 'pupil_mm');

    % Load optimal lambda from no-global results (reuse same lambda grid)
    f_res = fullfile(results_dir_nog, sprintf('subj%02d_results_no_global.mat', s));
    if ~exist(f_res, 'file')
        fprintf('  WARNING: no-global results not found for subject %d, skipping\n', s);
        continue;
    end
    load(f_res, 'Lambda', 'ind');

    % ----------------------------------------------------------------
    % Build X (intercept + global luminance only) and Y (pupil)
    % Same preprocessing as main analysis
    % ----------------------------------------------------------------
    ig = 120; d = 0;
    X_g = []; Y_g = [];

    for i = 1:length(global_screen)
        % X: intercept + global luminance column only
        X_g = [X_g; [ones(size(global_screen{i}(ig+1:end-d))) ...
                      global_screen{i}(ig+1:end-d)]];

        % Y: detrended, smoothed pupil
        y    = pupil_mm{i};
        nanx = isnan(y);
        t    = 1:numel(y);
        y(nanx) = interp1(t(~nanx), y(~nanx), t(nanx), 'linear', 'extrap');
        Yd   = detrend(y);
        ff   = normpdf(-10:10, 0, 0.1); ff = ff./sum(ff);
        Yd   = conv(Yd, ff, 'same');
        Y_g  = [Y_g; Yd((ig+1+d):end)];
    end

    % ----------------------------------------------------------------
    % Convert global luminance to pupil prediction via biomechanical model
    % Same pipeline as spatial model (so comparison is fair)
    % ----------------------------------------------------------------
    dt = 0.0167; taup = 0.15; taus = 0.64;
    P_g      = ones(size(X_g));   % column 1 = intercept (stays as ones)
    arrousal = zeros(size(X_g,1),1);
    SCinput  = zeros(size(X_g,1),1);

    t_vec = 1:numel(X_g(:,1));

    % Only column 2 (global luminance) — no patches
    nanx = isnan(X_g(:,2));
    X_g(nanx,2) = interp1(t_vec(~nanx), X_g(~nanx,2), t_vec(nanx), 'linear', 'extrap');
    [fs, fp]  = pupil_create_command(X_g(:,2)', arrousal', SCinput', taup, taus, dt);
    P_g(:,2) = pupil_model(fs, fp, dt);

    % ----------------------------------------------------------------
    % Fit regularized regression with cross-validated lambda
    % ----------------------------------------------------------------
    Lambda_grid = logspace(-5, 5, 30);

    CVMdl = fitrlinear(P_g, Y_g, 'ObservationsIn', 'rows', 'KFold', 5, ...
        'Lambda', Lambda_grid, 'Learner', 'leastsquares', ...
        'Solver', 'sparsa', 'Regularization', 'lasso');
    mse_g   = kfoldLoss(CVMdl);
    ind_g   = find(mse_g == min(mse_g), 1);

    Mdl_g = fitrlinear(P_g, Y_g, 'ObservationsIn', 'rows', ...
        'Lambda', Lambda_grid(ind_g), 'Learner', 'leastsquares', ...
        'Solver', 'sparsa', 'Regularization', 'lasso');

    % Compute R²
    Yp_g  = predict(Mdl_g, P_g);
    SSR_g = nansum((Yp_g - Y_g).^2);
    TSS_g = nansum((Y_g  - nanmean(Y_g)).^2);
    r2_global_only(idx) = 1 - SSR_g/TSS_g;

    fprintf('  Subject %02d: R² (global only) = %.4f\n', s, r2_global_only(idx));
end

% ====================================================================
% SUMMARY AND COMPARISON
% ====================================================================
fprintf('\n=============================================\n');
fprintf('RESULTS SUMMARY\n');
fprintf('=============================================\n');
fprintf('Global-only model:\n');
fprintf('  Mean R² = %.4f\n', nanmean(r2_global_only));
fprintf('  SD R²   = %.4f\n', nanstd(r2_global_only));
fprintf('  Range   = [%.4f, %.4f]\n', nanmin(r2_global_only), nanmax(r2_global_only));

% Load spatial-only R² for comparison
r2_spatial = nan(1, length(subjects));
for idx = 1:length(subjects)
    s = subjects(idx);
    f = fullfile(results_dir_nog, sprintf('subj%02d_results_no_global.mat', s));
    if exist(f, 'file')
        tmp = load(f, 'Rsquared');
        r2_spatial(idx) = tmp.Rsquared;
    end
end

valid = ~isnan(r2_global_only) & ~isnan(r2_spatial);
[~, p_comp, ~, stats_comp] = ttest(r2_spatial(valid)', r2_global_only(valid)');
mean_diff = mean(r2_spatial(valid) - r2_global_only(valid));
cohens_d  = mean_diff / std(r2_spatial(valid) - r2_global_only(valid));

fprintf('\nSpatial-only model:\n');
fprintf('  Mean R² = %.4f\n', nanmean(r2_spatial));
fprintf('  SD R²   = %.4f\n', nanstd(r2_spatial));

fprintf('\nComparison (spatial vs global-only):\n');
fprintf('  Mean delta-R²:  %.4f\n', mean_diff);
fprintf('  t(%d) = %.3f, p = %.4f\n', stats_comp.df, stats_comp.tstat, p_comp);
fprintf('  Cohen d = %.3f\n', cohens_d);

if mean_diff > 0 && p_comp < 0.05
    fprintf('\n  INTERPRETATION: Spatial model explains significantly MORE\n');
    fprintf('  variance than global luminance alone.\n');
    fprintf('  Local spatial structure contributes beyond global signal.\n');
elseif mean_diff > 0 && p_comp >= 0.05
    fprintf('\n  INTERPRETATION: Spatial model numerically outperforms\n');
    fprintf('  global-only but difference is not significant.\n');
else
    fprintf('\n  INTERPRETATION: No advantage for spatial model.\n');
end

% ====================================================================
% FIGURE: Paired scatter — spatial vs global-only R² per subject
% ====================================================================
fig = figure('Position', [100 100 500 460], 'Color', 'w');

scatter(r2_global_only(valid), r2_spatial(valid), 50, [0.2 0.5 0.8], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
hold on

ax_min = min([r2_global_only(valid) r2_spatial(valid)]) - 0.02;
ax_max = max([r2_global_only(valid) r2_spatial(valid)]) + 0.02;
plot([ax_min ax_max], [ax_min ax_max], 'k--', 'LineWidth', 1.5)

p_reg = polyfit(r2_global_only(valid), r2_spatial(valid), 1);
x_fit = linspace(ax_min, ax_max, 100);
plot(x_fit, polyval(p_reg, x_fit), 'r-', 'LineWidth', 2)

xlabel('R² (Global Luminance Only)',  'FontSize', 12, 'FontWeight', 'bold')
ylabel('R² (Spatial Patches Only)',   'FontSize', 12, 'FontWeight', 'bold')
title('Spatial vs Global-Only Model', 'FontSize', 13, 'FontWeight', 'bold')
axis equal
xlim([ax_min ax_max]); ylim([ax_min ax_max])
grid on

text(0.05, 0.97, sprintf('Mean ΔR² = %+.4f\nt(%d) = %.2f\np = %.4f\nd = %.3f', ...
    mean_diff, stats_comp.df, stats_comp.tstat, p_comp, cohens_d), ...
    'Units','normalized','FontSize',10,'VerticalAlignment','top','BackgroundColor','w','EdgeColor','k')

saveas(fig, fullfile(results_dir_out, 'FigS4_SpatialVsGlobalOnly.png'))
saveas(fig, fullfile(results_dir_out, 'FigS4_SpatialVsGlobalOnly.svg'))
close(fig)

% ====================================================================
% SAVE RESULTS
% ====================================================================
save(fullfile(results_dir_out, 'global_only_results.mat'), ...
    'r2_global_only', 'subjects', 'mean_diff', 'p_comp', 'stats_comp', 'cohens_d');

results_table = table(subjects(:), r2_global_only(:), r2_spatial(:), ...
    r2_spatial(:) - r2_global_only(:), ...
    'VariableNames', {'Subject','R2_GlobalOnly','R2_Spatial','Delta_R2'});
writetable(results_table, fullfile(results_dir_out, 'global_only_comparison.csv'));

fprintf('\n✓ Saved results to: %s\n', results_dir_out);
fprintf('  FigS4_SpatialVsGlobalOnly.png/.svg\n');
fprintf('  global_only_results.mat\n');
fprintf('  global_only_comparison.csv\n');
fprintf('=============================================\n');
