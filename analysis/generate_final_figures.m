% generate_final_figures.m
% Publication-ready figures - parula colormap throughouts

clear; close all;

% ====================================================================
% CONFIGURATION - UPDATE THESE PATHS
% ====================================================================
results_dir_noglobal = '/Users/jaspreetdodd/Desktop/results/results_no_global';
results_dir_global   = '/Users/jaspreetdodd/Desktop/results';
figures_dir          = fullfile(results_dir_noglobal, 'publication_figures_v2');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

good_subjects = [1:17 19:51]; % excluding a participant due to excessive data loss
subjects      = good_subjects;
% ====================================================================

fprintf('=============================================\n');
fprintf('GENERATING FINAL PUBLICATION FIGURES v2\n');
fprintf('N = %d subjects\n', length(subjects));
fprintf('=============================================\n\n');

%% ====================================================================
%  LOAD ALL DATA
%  ====================================================================
fprintf('Loading no-global results...\n');
all_weights  = [];
all_r2_nog   = [];
all_lambda   = [];
all_mse_nog  = [];
all_metadata = [];

for s = subjects
    f = fullfile(results_dir_noglobal, sprintf('subj%02d_results_no_global.mat', s));
    if exist(f, 'file')
        load(f, 'weight_map', 'Rsquared', 'Lambda', 'ind', 'mse', 'metadata');
        all_weights        = cat(3, all_weights, weight_map);
        all_r2_nog(end+1)  = Rsquared;
        all_lambda(end+1)  = Lambda(ind);
        all_mse_nog(end+1) = min(mse);
        if isempty(all_metadata), all_metadata = metadata; end
    else
        fprintf('  WARNING: Missing no-global data for subject %d\n', s);
    end
end
fprintf('  Loaded %d subjects (no-global).\n\n', size(all_weights,3));

fprintf('Loading with-global results...\n');
all_r2_glob = nan(1, length(subjects));
for idx = 1:length(subjects)
    s = subjects(idx);
    f = fullfile(results_dir_global, sprintf('subj%02d_results_compact.mat', s));
    if exist(f, 'file')
        load(f, 'Rsquared');
        all_r2_glob(idx) = Rsquared;
    else
        fprintf('  WARNING: Missing global data for subject %d\n', s);
    end
end
fprintf('  Loaded %d subjects (with-global).\n\n', sum(~isnan(all_r2_glob)));

nX    = all_metadata.nX;   % 62
nY    = all_metadata.nY;   % 49
nSubj = size(all_weights, 3);

%% ====================================================================
%  ECCENTRICITY MAP
%  ====================================================================
center_x     = nX / 2;
center_y     = nY / 2;
eccentricity = zeros(nY, nX);
for i = 1:nX
    for j = 1:nY
        eccentricity(j,i) = sqrt((i - center_x)^2 + (j - center_y)^2);
    end
end

%% ====================================================================
%  GROUP STATISTICS
%  ====================================================================
mean_weights = mean(all_weights, 3);
[~, p_patch] = ttest(all_weights, 0, 'dim', 3);

% 99th percentile colormap limit
clim_val = prctile(abs(all_weights(:)), 99);

% Central vs peripheral masks
central_mask    = eccentricity < 5;
peripheral_mask = eccentricity > 15 & eccentricity < 30;

central_w    = nan(nSubj,1);
peripheral_w = nan(nSubj,1);
for s = 1:nSubj
    w = all_weights(:,:,s);
    central_w(s)    = mean(w(central_mask),    'omitnan');
    peripheral_w(s) = mean(w(peripheral_mask), 'omitnan');
end

% Mean difference (central minus peripheral) - more interpretable than fold-change
mean_central    = mean(central_w);
mean_peripheral = mean(peripheral_w);
mean_diff_cp    = mean_central - mean_peripheral;

[~, p_paired, ~, stats_paired] = ttest(central_w, peripheral_w);
cohens_d = mean(central_w - peripheral_w) / std(central_w - peripheral_w);

% Vertical asymmetry
upper_mask = false(nY, nX);
lower_mask = false(nY, nX);
upper_mask(1:floor(nY/2), :)    = true;
lower_mask(ceil(nY/2)+1:end, :) = true;

upper_w = nan(nSubj,1);
lower_w = nan(nSubj,1);
for s = 1:nSubj
    w = all_weights(:,:,s);
    upper_w(s) = mean(w(upper_mask), 'omitnan');
    lower_w(s) = mean(w(lower_mask), 'omitnan');
end
[~, p_vert, ~, stats_vert] = ttest(upper_w, lower_w);
cohens_d_vert = mean(upper_w - lower_w) / std(upper_w - lower_w);

% Global vs no-global
valid_both     = ~isnan(all_r2_glob);
r2_nog_paired  = all_r2_nog(valid_both);
r2_glob_paired = all_r2_glob(valid_both);
[~, p_mc, ~, stats_mc] = ttest(r2_glob_paired, r2_nog_paired);
mean_diff_mc   = mean(r2_glob_paired - r2_nog_paired);
cohens_d_mc    = mean(r2_glob_paired - r2_nog_paired) / std(r2_glob_paired - r2_nog_paired);

fprintf('=============================================\n');
fprintf('KEY STATISTICS\n');
fprintf('=============================================\n');
fprintf('R² (no-global):   Mean=%.3f, SD=%.3f, Range=[%.3f %.3f]\n', ...
    mean(all_r2_nog), std(all_r2_nog), min(all_r2_nog), max(all_r2_nog));
fprintf('R² (with-global): Mean=%.3f, SD=%.3f\n', ...
    mean(r2_glob_paired), std(r2_glob_paired));
fprintf('Model comparison: ΔR²=%.4f, t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    mean_diff_mc, stats_mc.df, stats_mc.tstat, p_mc, cohens_d_mc);
fprintf('Central mean weight:    %.6f\n', mean_central);
fprintf('Peripheral mean weight: %.6f\n', mean_peripheral);
fprintf('Central-Peripheral diff: %.6f, p=%.2e, d=%.3f\n', mean_diff_cp, p_paired, cohens_d);
fprintf('Upper vs Lower: t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    stats_vert.df, stats_vert.tstat, p_vert, cohens_d_vert);
fprintf('=============================================\n\n');

%% ====================================================================
%  FIGURE 1: Representative Individual Weight Maps
%  ====================================================================
fprintf('Creating Figure 1: Representative weight maps...\n');

[~, idx_sorted] = sort(all_r2_nog, 'descend');
rep_idx = [
    idx_sorted(1);
    idx_sorted(round(length(idx_sorted)*0.25));
    idx_sorted(round(length(idx_sorted)*0.50));
    idx_sorted(round(length(idx_sorted)*0.75));
    idx_sorted(end-2);
    idx_sorted(end)
];

fig1 = figure('Position', [100 100 1800 900], 'Color', 'w');
for i = 1:6
    sidx = rep_idx(i);
    subplot(2,3,i)
    imagesc(all_weights(:,:,sidx))
    caxis([-clim_val clim_val])
    colormap(parula)
    if mod(i,3)==0
        c = colorbar;
        c.Label.String = 'Regression Weight';
        c.Label.FontSize = 10;
    end
    title(sprintf('Subject %d  (R² = %.3f)', subjects(sidx), all_r2_nog(sidx)), ...
        'FontSize', 11, 'FontWeight', 'bold')
    xlabel('Horizontal Patch (i)', 'FontSize', 9)
    ylabel('Vertical Patch (j)',   'FontSize', 9)
    axis equal tight
    hold on
    plot([nX/2 nX/2], [0.5 nY+0.5], 'w:', 'LineWidth', 1.2)
    plot([0.5 nX+0.5], [nY/2 nY/2], 'w:', 'LineWidth', 1.2)
end
sgtitle('Spatial Weight Maps: Representative Subjects', 'FontSize', 14, 'FontWeight', 'bold')

saveas(fig1, fullfile(figures_dir, 'Fig1_Representative_WeightMaps.png'))
saveas(fig1, fullfile(figures_dir, 'Fig1_Representative_WeightMaps.svg'))
close(fig1)
fprintf('  ✓ Saved Figure 1\n\n');

%% ====================================================================
%  FIGURE 2: Group Average + FDR-Corrected Significant Patches
%  ====================================================================
fprintf('Creating Figure 2: Group average + significant patches...\n');

% FDR correction
p_vec = p_patch(:);
[~, ~, ~, p_fdr_vec] = fdr_bh(p_vec, 0.05);
p_fdr = reshape(p_fdr_vec, nY, nX);

% Significant map: show mean weight where sig, NaN elsewhere
sig_map = mean_weights;
sig_map(p_fdr > 0.05) = NaN;
n_sig = sum(~isnan(sig_map(:)));

% Build RGB image: grey for non-significant, parula for significant
cmap     = parula(256);
grey_rgb = [0.78 0.78 0.78];

sig_norm = (sig_map - (-clim_val)) / (2 * clim_val);
sig_norm = min(max(sig_norm, 0), 1);

rgb_img = repmat(reshape(grey_rgb, [1 1 3]), nY, nX);
for row = 1:nY
    for col = 1:nX
        if ~isnan(sig_map(row,col))
            cidx = max(1, min(256, round(sig_norm(row,col) * 255) + 1));
            rgb_img(row,col,:) = cmap(cidx,:);
        end
    end
end

fig2 = figure('Position', [100 100 1400 520], 'Color', 'w');

% Panel A: Group average
subplot(1,2,1)
imagesc(mean_weights)
caxis([-clim_val clim_val])
colormap(gca, parula)
c = colorbar;
c.Label.String = 'Mean Regression Weight';
c.Label.FontSize = 11;
title(sprintf('A. Group Average  (N=%d)', nSubj), 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Horizontal Patch (i)', 'FontSize', 11)
ylabel('Vertical Patch (j)',   'FontSize', 11)
axis equal tight
hold on
plot([nX/2 nX/2], [0.5 nY+0.5], 'w:', 'LineWidth', 1.5)
plot([0.5 nX+0.5], [nY/2 nY/2], 'w:', 'LineWidth', 1.5)

% Panel B: FDR-corrected significant patches with grey masking
subplot(1,2,2)
image(rgb_img)
axis equal tight
hold on
plot([nX/2 nX/2], [0.5 nY+0.5], 'w:', 'LineWidth', 1.5)
plot([0.5 nX+0.5], [nY/2 nY/2], 'w:', 'LineWidth', 1.5)
colormap(gca, parula)
clim([-clim_val clim_val])
c2 = colorbar;
c2.Label.String = 'Mean Weight (FDR q<0.05)';
c2.Label.FontSize = 11;

% Grey patch legend item
annotation('rectangle', [0.895 0.14 0.018 0.04], ...
    'FaceColor', grey_rgb, 'EdgeColor', [0.4 0.4 0.4])
annotation('textbox', [0.916 0.125 0.06 0.06], ...
    'String', 'n.s.', 'FontSize', 9, 'EdgeColor', 'none')

title(sprintf('B. FDR-Corrected Significant Patches\n(%d / %d, q<0.05)', n_sig, nX*nY), ...
    'FontSize', 13, 'FontWeight', 'bold')
xlabel('Horizontal Patch (i)', 'FontSize', 11)
ylabel('Vertical Patch (j)',   'FontSize', 11)

saveas(fig2, fullfile(figures_dir, 'Fig2_Group_Average.png'))
saveas(fig2, fullfile(figures_dir, 'Fig2_Group_Average.svg'))
close(fig2)
fprintf('  ✓ Saved Figure 2  (%d/%d patches significant after FDR)\n\n', n_sig, nX*nY);

%% ====================================================================
%  FIGURE 3: Model Performance
%  ====================================================================
fprintf('Creating Figure 3: Model performance...\n');

fig3 = figure('Position', [100 100 1800 460], 'Color', 'w');

% Panel A: R² distribution
subplot(1,4,1)
histogram(all_r2_nog, 15, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.8)
hold on
yl = ylim;
plot([mean(all_r2_nog)   mean(all_r2_nog)],   yl, 'r--', 'LineWidth', 2)
plot([median(all_r2_nog) median(all_r2_nog)], yl, 'g--', 'LineWidth', 2)
xlabel('R²', 'FontSize', 11)
ylabel('Number of Subjects', 'FontSize', 11)
title(sprintf('A. R² Distribution\n(N=%d)', nSubj), 'FontSize', 12, 'FontWeight', 'bold')
legend({'Data','Mean','Median'}, 'Location','northwest', 'FontSize', 9)
grid on
text(0.97, 0.97, sprintf('μ = %.3f\nσ = %.3f\nMdn = %.3f', ...
    mean(all_r2_nog), std(all_r2_nog), median(all_r2_nog)), ...
    'Units','normalized','FontSize',9,'VerticalAlignment','top', ...
    'HorizontalAlignment','right','BackgroundColor','w')

% Panel B: Individual R² sorted
subplot(1,4,2)
sorted_r2 = sort(all_r2_nog, 'descend');
bar(1:length(sorted_r2), sorted_r2, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none')
hold on
plot([0 length(sorted_r2)+1], [mean(all_r2_nog) mean(all_r2_nog)], 'r--', 'LineWidth', 2)
xlabel('Subject (sorted by R²)', 'FontSize', 11)
ylabel('R²', 'FontSize', 11)
title('B. Individual Performance', 'FontSize', 12, 'FontWeight', 'bold')
grid on
ylim([0 max(sorted_r2)*1.1])
xlim([0 length(sorted_r2)+1])

% Panel C: Central vs Peripheral
% Report absolute means and mean difference, not fold-change
subplot(1,4,3)
boxplot([central_w, peripheral_w], ...
    'Labels', {'Central (0-5°)', 'Peripheral (15-30°)'}, ...
    'Widths', 0.5)
ylabel('Mean Regression Weight', 'FontSize', 11)
title('C. Central vs. Peripheral', 'FontSize', 12, 'FontWeight', 'bold')
grid on
hold on
% Individual subject lines (light grey)
for s = 1:nSubj
    plot([1 2], [central_w(s) peripheral_w(s)], ...
        'Color', [0.6 0.6 0.6 0.2], 'LineWidth', 0.5)
end
% Stats annotation — mean difference instead of fold-change
text(0.5, 0.97, sprintf('Δ = %.4f\np = %.2e\nd = %.2f', ...
    mean_diff_cp, p_paired, cohens_d), ...
    'Units','normalized','FontSize',9,'VerticalAlignment','top', ...
    'HorizontalAlignment','center','BackgroundColor','w')

% Panel D: Global vs No-Global scatter
subplot(1,4,4)
scatter(r2_nog_paired, r2_glob_paired, 45, [0.2 0.5 0.8], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
hold on
ax_min = min([r2_nog_paired r2_glob_paired]) - 0.02;
ax_max = max([r2_nog_paired r2_glob_paired]) + 0.02;
plot([ax_min ax_max], [ax_min ax_max], 'k--', 'LineWidth', 1.5)  % unity line
p_reg = polyfit(r2_nog_paired, r2_glob_paired, 1);
x_fit = linspace(ax_min, ax_max, 100);
plot(x_fit, polyval(p_reg, x_fit), 'r-', 'LineWidth', 2)
xlabel('R² (Spatial Only)',             'FontSize', 11)
ylabel('R² (Spatial + Global Lum.)',    'FontSize', 11)
title('D. Model Comparison', 'FontSize', 12, 'FontWeight', 'bold')
axis equal
xlim([ax_min ax_max]); ylim([ax_min ax_max])
grid on
text(0.05, 0.97, sprintf('Mean ΔR² = %+.4f\nt(%d) = %.2f\np = %.4f\nd = %.3f', ...
    mean_diff_mc, stats_mc.df, stats_mc.tstat, p_mc, cohens_d_mc), ...
    'Units','normalized','FontSize',9,'VerticalAlignment','top','BackgroundColor','w')

saveas(fig3, fullfile(figures_dir, 'Fig3_Model_Performance.png'))
saveas(fig3, fullfile(figures_dir, 'Fig3_Model_Performance.svg'))
close(fig3)
fprintf('  ✓ Saved Figure 3\n\n');

%% ================================
%  FIGURE 4: Radial Weight Profile 
%  =================================
fprintf('Creating Figure 4: Radial weight profile...\n');

bins        = [0 2 5 10 15 20 30];
bin_centers = (bins(1:end-1) + bins(2:end)) / 2;
nBins       = numel(bins) - 1;

mean_by_ecc = nan(nBins, nSubj);
for s = 1:nSubj
    w = all_weights(:,:,s);
    for b = 1:nBins
        mask = eccentricity >= bins(b) & eccentricity < bins(b+1);
        mean_by_ecc(b,s) = mean(w(mask), 'omitnan');
    end
end

grand_mean = mean(mean_by_ecc, 2, 'omitnan');
grand_se   = std(mean_by_ecc,  0, 2, 'omitnan') / sqrt(nSubj);
grand_ci   = 1.96 * grand_se;

% Per-bin t-tests with FDR
p_bin = nan(nBins,1);
for b = 1:nBins
    y = mean_by_ecc(b,:);
    y = y(~isnan(y));
    if numel(y) >= 3
        [~, p_bin(b)] = ttest(y, 0);
    end
end
[~, ~, ~, p_bin_fdr] = fdr_bh(p_bin, 0.05);

% Subject-level exponential fits for INFERENCE
ft     = fittype('a*exp(-b*x) + c');
b_subj = nan(nSubj,1);
for s = 1:nSubj
    y  = mean_by_ecc(:,s);
    ok = ~isnan(y);
    if sum(ok) >= 4
        try
            sp   = [max(y(ok)), 0.1, min(y(ok))];
            fsub = fit(bin_centers(ok)', y(ok), ft, 'StartPoint', sp, ...
                'Lower', [-Inf, 0, -Inf]);
            b_subj(s) = fsub.b;
        catch
            b_subj(s) = NaN;
        end
    end
end
b_mean = mean(b_subj, 'omitnan');
b_se   = std(b_subj,  'omitnan') / sqrt(sum(~isnan(b_subj)));
[~, p_decay, ~, stats_decay] = ttest(b_subj, 0, 'Tail', 'right');

% Quadratic fit for VISUALIZATION
x_fit  = linspace(0, max(bins), 200);
p_quad = polyfit(bin_centers', grand_mean, 2);
y_quad = polyval(p_quad, x_fit);
y_pred_bins = polyval(p_quad, bin_centers');
r2_quad = 1 - sum((grand_mean - y_pred_bins).^2) / ...
               sum((grand_mean - mean(grand_mean)).^2);

fig4 = figure('Position', [100 100 800 620], 'Color', 'w');

% Main plot — leave extra room at top for stars by setting ylim explicitly
ax = axes;
errorbar(bin_centers, grand_mean, grand_ci, 'o-', ...
    'LineWidth', 2.5, 'MarkerSize', 9, 'CapSize', 10, ...
    'Color', [0.15 0.45 0.75], 'MarkerFaceColor', [0.15 0.45 0.75])
hold on
plot([0 max(bins)], [0 0], 'k--', 'LineWidth', 1.5)
plot(x_fit, y_quad, 'r--', 'LineWidth', 2)

% Set ylim first so we can place stars cleanly above data
data_top = max(grand_mean + grand_ci);
ylim([min(grand_mean - grand_ci) * 1.3, data_top * 1.35])

% Significance stars — placed at fixed y above the top of error bars
y_star = data_top * 1.18;
for b = 1:nBins
    if ~isnan(p_bin_fdr(b)) && p_bin_fdr(b) < 0.05
        text(bin_centers(b), y_star, '*', ...
            'HorizontalAlignment','center','FontSize',18,'FontWeight','bold', ...
            'Color', [0.1 0.1 0.1])
    end
end

xlabel('Eccentricity from Fixation (°)', 'FontSize', 13, 'FontWeight', 'bold')
ylabel('Mean Regression Weight',          'FontSize', 13, 'FontWeight', 'bold')
title('Pupil Sensitivity by Retinal Eccentricity', ...
    'FontSize', 14, 'FontWeight', 'bold')   % no N= in title (put in legend)
grid on
xlim([-0.5 max(bins)+0.5])

legend({sprintf('Mean ± 95% CI  (N=%d)', nSubj), ...
        'Zero baseline', ...
        sprintf('Quadratic fit (R²=%.3f)', r2_quad)}, ...
    'Location','northeast','FontSize',10)

% Stats box — inferential results from subject-level exponential fits
text(0.97, 0.62, sprintf(['Peak weight: %.4f\n' ...
                           'Exp decay b: %.3f ± %.3f /°\n' ...
                           'Test b>0: p = %.2e\n' ...
                           'Central mean: %.4f\n' ...
                           'Peripheral mean: %.4f'], ...
    grand_mean(1), b_mean, b_se, p_decay, mean_central, mean_peripheral), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',10, ...
    'BackgroundColor','w','EdgeColor','k','Margin',5)

saveas(fig4, fullfile(figures_dir, 'Fig4_Radial_Profile.png'))
saveas(fig4, fullfile(figures_dir, 'Fig4_Radial_Profile.svg'))
close(fig4)
fprintf('  ✓ Saved Figure 4\n\n');

%% ====================================================================
%  FIGURE 5: Reliability Analysis — FIXED consistency colormap
%  ====================================================================
fprintf('Creating Figure 5: Reliability analysis...\n');

fig5 = figure('Position', [100 100 1400 480], 'Color', 'w');

% Panel A: Split-half reliability
subplot(1,3,1)
odd_idx  = 1:2:nSubj;
even_idx = 2:2:nSubj;
mean_odd  = mean(all_weights(:,:,odd_idx),  3);
mean_even = mean(all_weights(:,:,even_idx), 3);

scatter(mean_odd(:), mean_even(:), 4, [0.3 0.6 0.9], ...
    'filled', 'MarkerFaceAlpha', 0.35)
hold on
p_sh = polyfit(mean_odd(:), mean_even(:), 1);
x_sh = linspace(min(mean_odd(:)), max(mean_odd(:)), 100);
plot(x_sh, polyval(p_sh, x_sh), 'r-', 'LineWidth', 2)
plot(x_sh, x_sh, 'k--', 'LineWidth', 1)
[r_split, ~] = corr(mean_odd(:), mean_even(:));
xlabel('Odd Subjects Mean Weight',  'FontSize', 11)
ylabel('Even Subjects Mean Weight', 'FontSize', 11)
title('A. Split-Half Reliability',  'FontSize', 12, 'FontWeight', 'bold')
axis equal; grid on
text(0.05, 0.95, sprintf('r = %.3f\np < 0.001', r_split), ...
    'Units','normalized','FontSize',10,'VerticalAlignment','top','BackgroundColor','w')

% Panel B: Spatial consistency map — FIXED
% High consistency (1.0) = bright yellow, low/chance (0.5) = dark blue
% Use parula with floor at 0.5 so the central high-consistency region is warm
subplot(1,3,2)
consistency_map = sum(sign(all_weights) == sign(mean_weights), 3) / nSubj;

imagesc(consistency_map, [0.5 1])
colormap(gca, parula)   % parula: dark blue at low end, yellow at high end
c = colorbar;
c.Label.String = 'Proportion Same Sign as Group';
c.Label.FontSize = 10;
c.Ticks      = [0.5 0.6 0.7 0.8 0.9 1.0];
c.TickLabels = {'0.5 (chance)', '0.6', '0.7', '0.8', '0.9', '1.0'};

title('B. Spatial Consistency',   'FontSize', 12, 'FontWeight', 'bold')
xlabel('Horizontal Patch (i)',    'FontSize', 11)
ylabel('Vertical Patch (j)',      'FontSize', 11)
axis equal tight
hold on
plot([nX/2 nX/2], [0.5 nY+0.5], 'w:', 'LineWidth', 1.5)
plot([0.5 nX+0.5], [nY/2 nY/2], 'w:', 'LineWidth', 1.5)

% Panel C: Bootstrap reliability
subplot(1,3,3)
n_boot = 1000;
boot_r = zeros(n_boot,1);
fprintf('  Running bootstrap (N=%d)...\n', n_boot);
for b = 1:n_boot
    bidx      = randi(nSubj, [nSubj 1]);
    boot_mean = mean(all_weights(:,:,bidx), 3);
    boot_r(b) = corr(boot_mean(:), mean_weights(:));
    if mod(b,250)==0, fprintf('    %d/%d\n', b, n_boot); end
end
ci_lo = prctile(boot_r, 2.5);
ci_hi = prctile(boot_r, 97.5);

histogram(boot_r, 30, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k', 'FaceAlpha', 0.8)
xlabel('Correlation with Full Sample', 'FontSize', 11)
ylabel('Frequency',                    'FontSize', 11)
title('C. Bootstrap Reliability',      'FontSize', 12, 'FontWeight', 'bold')
grid on
text(0.05, 0.95, sprintf('Mean r = %.3f\n95%% CI: [%.3f, %.3f]', ...
    mean(boot_r), ci_lo, ci_hi), ...
    'Units','normalized','FontSize',9,'VerticalAlignment','top','BackgroundColor','w')

saveas(fig5, fullfile(figures_dir, 'Fig5_Reliability.png'))
saveas(fig5, fullfile(figures_dir, 'Fig5_Reliability.svg'))
close(fig5)
fprintf('  ✓ Saved Figure 5\n\n');

%% ====================================================================
%  STATISTICS TABLE
%  ====================================================================
fprintf('Saving statistics table...\n');

stats_out        = table();
stats_out.Metric = {
    'N subjects analyzed';
    'N excluded';
    'Mean R² (spatial only)';
    'SD R² (spatial only)';
    'Median R² (spatial only)';
    'Min R²';
    'Max R²';
    'Mean R² (spatial + global)';
    'SD R² (spatial + global)';
    'Mean delta-R² (global adds)';
    'Model comparison t-stat (df)';
    'Model comparison p-value';
    'Model comparison Cohen d';
    'Central mean weight';
    'Peripheral mean weight';
    'Central-Peripheral mean difference';
    'Central vs Peripheral p-value';
    'Central vs Peripheral Cohen d';
    'Upper vs Lower VF t-stat (df)';
    'Upper vs Lower VF p-value';
    'Upper vs Lower VF Cohen d';
    'Exp decay b mean';
    'Exp decay b SE';
    'Exp decay b>0 p-value';
    'N patches significant (FDR q<0.05)';
    'Total patches';
    'Proportion significant';
    'Split-half r';
    'Bootstrap mean r';
    'Bootstrap 95% CI low';
    'Bootstrap 95% CI high';
    'Optimal lambda mean';
    'Optimal lambda median'
};
stats_out.Value = {
    sprintf('%d', nSubj);
    '1 (subject 18, >20% missing data)';
    sprintf('%.4f', mean(all_r2_nog));
    sprintf('%.4f', std(all_r2_nog));
    sprintf('%.4f', median(all_r2_nog));
    sprintf('%.4f', min(all_r2_nog));
    sprintf('%.4f', max(all_r2_nog));
    sprintf('%.4f', mean(r2_glob_paired));
    sprintf('%.4f', std(r2_glob_paired));
    sprintf('%.4f', mean_diff_mc);
    sprintf('%.3f (df=%d)', stats_mc.tstat, stats_mc.df);
    sprintf('%.4f', p_mc);
    sprintf('%.3f', cohens_d_mc);
    sprintf('%.6f', mean_central);
    sprintf('%.6f', mean_peripheral);
    sprintf('%.6f', mean_diff_cp);
    sprintf('%.2e', p_paired);
    sprintf('%.3f', cohens_d);
    sprintf('%.3f (df=%d)', stats_vert.tstat, stats_vert.df);
    sprintf('%.4f', p_vert);
    sprintf('%.3f', cohens_d_vert);
    sprintf('%.4f', b_mean);
    sprintf('%.4f', b_se);
    sprintf('%.2e', p_decay);
    sprintf('%d', n_sig);
    sprintf('%d', nX*nY);
    sprintf('%.3f', n_sig/(nX*nY));
    sprintf('%.3f', r_split);
    sprintf('%.3f', mean(boot_r));
    sprintf('%.3f', ci_lo);
    sprintf('%.3f', ci_hi);
    sprintf('%.6f', mean(all_lambda));
    sprintf('%.6f', median(all_lambda))
};

writetable(stats_out, fullfile(figures_dir, 'Statistics_Complete.csv'));
fprintf('  ✓ Saved Statistics_Complete.csv\n\n');

%% ====================================================================
%  DONE
%  ====================================================================
fprintf('=============================================\n');
fprintf('ALL FIGURES GENERATED SUCCESSFULLY\n');
fprintf('Saved to: %s\n', figures_dir);
fprintf('=============================================\n');


%% ====================================================================
%  Now run additional analyses
%  ====================================================================
% Run AFTER above (relies on workspace variables)
% Adds: effect size for significant patches, individual differences,
%       permutation test, and vertical asymmetry figure annotation.
%
% Requires in workspace: all_weights, all_r2_nog, mean_weights,
%   p_fdr, sig_map, n_sig, nSubj, nX, nY, subjects, figures_dir,
%   results_dir_noglobal, upper_w, lower_w, p_vert, stats_vert,
%   cohens_d_vert, eccentricity, all_metadata

fprintf('\n=============================================\n');
fprintf('ADDITIONAL ANALYSES\n');
fprintf('=============================================\n\n');
%% ====================================================================
%  ITEM 2: Effect size for FDR-significant patches
%  ====================================================================
fprintf('Item 2: Effect size for significant patches...\n');

% For each significant patch, compute Cohen d across subjects
% (mean weight / SD of weights across subjects)
patch_mean = mean(all_weights, 3);
patch_sd   = std(all_weights,  0, 3);
patch_d    = patch_mean ./ patch_sd;   % Cohen d per patch

% Restrict to FDR-significant patches
sig_d = patch_d(~isnan(sig_map));   % sig_map has NaN for non-sig patches

mean_d_sig    = mean(abs(sig_d));
median_d_sig  = median(abs(sig_d));
mean_w_sig    = mean(abs(patch_mean(~isnan(sig_map))));

% Also: proportion of significant patches that are POSITIVE vs NEGATIVE
n_sig_pos = sum(patch_mean(~isnan(sig_map)) > 0);
n_sig_neg = sum(patch_mean(~isnan(sig_map)) < 0);

fprintf('  Significant patches: %d total\n', n_sig);
fprintf('  Positive: %d  |  Negative: %d\n', n_sig_pos, n_sig_neg);
fprintf('  Mean |Cohen d| across sig patches: %.3f\n', mean_d_sig);
fprintf('  Median |Cohen d|:                  %.3f\n', median_d_sig);
fprintf('  Mean |weight| in sig patches:      %.6f\n', mean_w_sig);

%% ====================================================================
%  ITEM 3: Model comparison interpretation — print clean summary
%  ====================================================================
fprintf('\nItem 3: Model comparison summary...\n');
% (stats already computed in main script)
% Just print a clean interpretation-ready summary
fprintf('  Spatial-only model:         R² = %.4f ± %.4f\n', ...
    mean(all_r2_nog), std(all_r2_nog));
fprintf('  Spatial + global model:     R² = %.4f ± %.4f\n', ...
    mean(r2_glob_paired), std(r2_glob_paired));
fprintf('  Mean delta-R²:              %.4f\n', mean_diff_mc);
fprintf('  Paired t-test:              t(%d) = %.3f, p = %.4f\n', ...
    stats_mc.df, stats_mc.tstat, p_mc);
fprintf('  Cohen d:                    %.3f\n', cohens_d_mc);
if p_mc > 0.05
    fprintf('  INTERPRETATION: Adding global luminance does NOT\n');
    fprintf('  significantly improve prediction. Spatial patches\n');
    fprintf('  subsume the global signal.\n');
else
    fprintf('  INTERPRETATION: Adding global luminance significantly\n');
    fprintf('  improves prediction.\n');
end

%% ====================================================================
%  ITEM 4: Individual differences — R² vs proportion missing data
%  ====================================================================
fprintf('\nItem 4: Individual differences (R² vs missing data)...\n');

prop_missing = nan(1, nSubj);

for idx = 1:nSubj
    s = subjects(idx);
    f = fullfile(results_dir_noglobal, sprintf('subj%02d_results_no_global.mat', s));
    if exist(f, 'file')
        % Load the model object to get access to the data indirectly
        % We compute missing rate from the retinalcoor file
        f2 = fullfile('/Users/jaspreetdodd/Desktop/local analysis/data', ...
            sprintf('subj%02d_retinalcoor.mat', s));
        if exist(f2, 'file')
            tmp = load(f2, 'pupil_mm');
            all_pupil = cat(1, tmp.pupil_mm{:});
            prop_missing(idx) = mean(isnan(all_pupil));
        end
    end
end

% Correlation between R² and proportion missing
valid_idx = ~isnan(prop_missing);
[r_miss, p_miss] = corr(all_r2_nog(valid_idx)', prop_missing(valid_idx)', ...
    'Type', 'Spearman');

fprintf('  Subjects with missing data computed: %d\n', sum(valid_idx));
fprintf('  Spearman r (R² vs prop missing): %.3f, p = %.4f\n', r_miss, p_miss);

% Figure: scatter of R² vs proportion missing
fig_id = figure('Position', [100 100 500 420], 'Color', 'w');
scatter(prop_missing(valid_idx)*100, all_r2_nog(valid_idx), ...
    50, [0.2 0.5 0.8], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
hold on
% Regression line
p_reg_miss = polyfit(prop_missing(valid_idx)', all_r2_nog(valid_idx)', 1);
x_miss     = linspace(min(prop_missing(valid_idx)), max(prop_missing(valid_idx)), 100);
plot(x_miss*100, polyval(p_reg_miss, x_miss), 'r-', 'LineWidth', 2)

xlabel('Missing Data (%)',  'FontSize', 12, 'FontWeight', 'bold')
ylabel('R²',                'FontSize', 12, 'FontWeight', 'bold')
title('Individual Differences: Model Fit vs Data Quality', ...
    'FontSize', 12, 'FontWeight', 'bold')
grid on
text(0.97, 0.97, sprintf('Spearman r = %.3f\np = %.4f', r_miss, p_miss), ...
    'Units','normalized','FontSize',11,'VerticalAlignment','top', ...
    'HorizontalAlignment','right','BackgroundColor','w','EdgeColor','k')

saveas(fig_id, fullfile(figures_dir, 'FigS1_IndividualDifferences.png'))
saveas(fig_id, fullfile(figures_dir, 'FigS1_IndividualDifferences.svg'))
close(fig_id)
fprintf('  ✓ Saved FigS1_IndividualDifferences\n');

%% ====================================================================
%  ITEM 5: Vertical asymmetry — annotate and save dedicated figure
%  ====================================================================
fprintf('\nItem 5: Vertical asymmetry...\n');

fprintf('  Upper VF mean weight: %.6f\n', mean(upper_w));
fprintf('  Lower VF mean weight: %.6f\n', mean(lower_w));
fprintf('  t(%d) = %.3f, p = %.4f, d = %.3f\n', ...
    stats_vert.df, stats_vert.tstat, p_vert, cohens_d_vert);

if mean(upper_w) > mean(lower_w)
    fprintf('  DIRECTION: Upper VF > Lower VF (consistent with literature)\n');
else
    fprintf('  DIRECTION: Lower VF > Upper VF\n');
end

% Simple bar plot for VF asymmetry
fig_va = figure('Position', [100 100 400 420], 'Color', 'w');
data_va = [upper_w, lower_w];
% Individual subject lines first
for s = 1:nSubj
    plot([1 2], [upper_w(s) lower_w(s)], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.8)
    hold on
end
boxplot(data_va, 'Labels', {'Upper VF', 'Lower VF'}, 'Widths', 0.4)
ylabel('Mean Regression Weight', 'FontSize', 12, 'FontWeight', 'bold')
title('Vertical Visual Field Asymmetry', 'FontSize', 12, 'FontWeight', 'bold')
grid on

% Significance bracket
y_top = max([upper_w; lower_w]) * 1.15;
plot([1 2], [y_top y_top], 'k-', 'LineWidth', 1.2)
if p_vert < 0.001
    sig_str = '***';
elseif p_vert < 0.01
    sig_str = '**';
elseif p_vert < 0.05
    sig_str = '*';
else
    sig_str = sprintf('p = %.3f', p_vert);
end
text(1.5, y_top * 1.02, sig_str, 'HorizontalAlignment','center','FontSize',13)

text(0.5, 0.05, sprintf('t(%d) = %.2f\np = %.4f\nd = %.3f', ...
    stats_vert.df, stats_vert.tstat, p_vert, cohens_d_vert), ...
    'Units','normalized','FontSize',10,'BackgroundColor','w','EdgeColor','k')

saveas(fig_va, fullfile(figures_dir, 'FigS2_VerticalAsymmetry.png'))
saveas(fig_va, fullfile(figures_dir, 'FigS2_VerticalAsymmetry.svg'))
close(fig_va)
fprintf('  ✓ Saved FigS2_VerticalAsymmetry\n');

%% ====================================================================
%  ITEM 6: Permutation test (all subjects)
%  ====================================================================
fprintf('\nItem 6: Permutation test...\n');
fprintf('  Running permutation for all %d subjects...\n', nSubj);
addpath('/Users/jaspreetdodd/Desktop/local analysis/code/actually used');

n_perms     = 100;   % 100 permutations per subject is sufficient for sanity check
perm_r2     = nan(nSubj, n_perms);

% We need access to the data files for this
% Load each subject's data, shuffle pupil, refit model
for idx = 1:nSubj
    s  = subjects(idx);
    f2 = fullfile('/Users/jaspreetdodd/Desktop/local analysis/data', ...
        sprintf('subj%02d_retinalcoor.mat', s));
    f_res = fullfile(results_dir_noglobal, ...
        sprintf('subj%02d_results_no_global.mat', s));

    if ~exist(f2, 'file') || ~exist(f_res, 'file'), continue; end

    % Load retinalcoor data
    tmp = load(f2, 'matrix', 'global_screen', 'pupil_mm', 'metadata');

    % Reconstruct X and Y exactly as in main analysis
    ig = 120; d = 0;
    X_s = []; Y_s = [];
    for i = 1:length(tmp.matrix)
        X_s = [X_s; [ones(size(tmp.global_screen{i}(ig+1:end-d))) ...
                      tmp.matrix{i}(ig+1:end-d,:)]];
        y   = tmp.pupil_mm{i};
        nanx = isnan(y); t_idx = 1:numel(y);
        y(nanx) = interp1(t_idx(~nanx), y(~nanx), t_idx(nanx), 'linear', 'extrap');
        Yd = detrend(y);
        ff = normpdf(-10:10, 0, 0.1); ff = ff./sum(ff);
        Yd = conv(Yd, ff, 'same');
        Y_s = [Y_s; Yd((ig+1+d):end)];
    end

    % Convert luminance to pupil prediction (same params as main analysis)
    dt = 0.0167; taup = 0.15; taus = 0.64;
    P_s = ones(size(X_s));
    arrousal = zeros(size(X_s,1),1);
    SCinput  = zeros(size(X_s,1),1);
    t_vec = 1:numel(X_s(:,1));
    for col = 2:size(X_s,2)
        nanx = isnan(X_s(:,col));
        X_s(nanx,col) = interp1(t_vec(~nanx), X_s(~nanx,col), ...
            t_vec(nanx), 'linear', 'extrap');
        [fs, fp]  = pupil_create_command(X_s(:,col)', arrousal', SCinput', taup, taus, dt);
        P_s(:,col) = pupil_model(fs, fp, dt);
    end

    % Load optimal lambda from saved results
    res = load(f_res, 'Lambda', 'ind');
    opt_lambda = res.Lambda(res.ind);

    % Run permutations
    for perm = 1:n_perms
        Y_perm = Y_s(randperm(length(Y_s)));   % shuffle pupil
        Mdl_perm = fitrlinear(P_s, Y_perm, 'ObservationsIn', 'rows', ...
            'Lambda', opt_lambda, 'Learner', 'leastsquares', ...
            'Solver', 'sparsa', 'Regularization', 'lasso');
        Yp_perm = predict(Mdl_perm, P_s);
        SSR_p   = nansum((Yp_perm - Y_perm).^2);
        TSS_p   = nansum((Y_perm  - nanmean(Y_perm)).^2);
        perm_r2(idx, perm) = 1 - SSR_p/TSS_p;
    end

    if mod(idx, 10) == 0
        fprintf('    Completed %d/%d subjects\n', idx, nSubj);
    end
end

% Summary
mean_perm_r2 = nanmean(perm_r2(:));
sd_perm_r2   = nanstd(perm_r2(:));
mean_subj_perm = nanmean(perm_r2, 2);   % per-subject mean permuted R²

fprintf('  Permutation R² (mean ± SD): %.4f ± %.4f\n', mean_perm_r2, sd_perm_r2);
fprintf('  True R²       (mean ± SD):  %.4f ± %.4f\n', mean(all_r2_nog), std(all_r2_nog));
fprintf('  True R² is %.1fx higher than permuted R²\n', mean(all_r2_nog)/mean_perm_r2);

% Figure: true vs permuted R² per subject
fig_perm = figure('Position', [100 100 600 440], 'Color', 'w');
hold on
% Permuted distribution as box
boxplot(mean_subj_perm, 'Positions', 1, 'Widths', 0.4, 'Colors', [0.7 0.7 0.7])
% True R² as scatter
scatter(ones(nSubj,1)*2, all_r2_nog', 50, [0.2 0.5 0.8], ...
    'filled', 'MarkerEdgeColor','k','LineWidth',0.5)
plot([1.8 2.2], [mean(all_r2_nog) mean(all_r2_nog)], 'r-', 'LineWidth', 2.5)

xlim([0.5 2.5])
ylim([min([mean_subj_perm; all_r2_nog'])*0.8, max(all_r2_nog)*1.15])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Permuted', 'True'}, 'FontSize', 11)
ylabel('R²', 'FontSize', 12, 'FontWeight', 'bold')
title('Permutation Sanity Check', 'FontSize', 12, 'FontWeight', 'bold')
grid on

text(0.97, 0.97, sprintf('Permuted: %.4f ± %.4f\nTrue: %.4f ± %.4f', ...
    mean_perm_r2, sd_perm_r2, mean(all_r2_nog), std(all_r2_nog)), ...
    'Units','normalized','FontSize',10,'VerticalAlignment','top', ...
    'HorizontalAlignment','right','BackgroundColor','w','EdgeColor','k')

saveas(fig_perm, fullfile(figures_dir, 'FigS3_PermutationTest.png'))
saveas(fig_perm, fullfile(figures_dir, 'FigS3_PermutationTest.svg'))
close(fig_perm)
fprintf('  ✓ Saved FigS3_PermutationTest\n');

%% ====================================================================
%  UPDATED STATISTICS TABLE — append new stats
%  ====================================================================
fprintf('\nUpdating statistics table...\n');

new_rows = {
    'Mean |Cohen d| in sig patches',     sprintf('%.3f', mean_d_sig);
    'Median |Cohen d| in sig patches',   sprintf('%.3f', median_d_sig);
    'Mean |weight| in sig patches',      sprintf('%.6f', mean_w_sig);
    'N significant patches (positive)',  sprintf('%d', n_sig_pos);
    'N significant patches (negative)',  sprintf('%d', n_sig_neg);
    'Spearman r (R² vs missing data)',   sprintf('%.3f', r_miss);
    'p-value (R² vs missing data)',      sprintf('%.4f', p_miss);
    'Upper VF mean weight',              sprintf('%.6f', mean(upper_w));
    'Lower VF mean weight',              sprintf('%.6f', mean(lower_w));
    'Permuted R² mean',                  sprintf('%.4f', mean_perm_r2);
    'Permuted R² SD',                    sprintf('%.4f', sd_perm_r2);
    'True/Permuted R² ratio',            sprintf('%.1fx', mean(all_r2_nog)/mean_perm_r2);
};

new_table = table(new_rows(:,1), new_rows(:,2), ...
    'VariableNames', {'Metric','Value'});

% Load existing table and append
existing = readtable(fullfile(figures_dir, 'Statistics_Complete.csv'));
combined = [existing; new_table];
writetable(combined, fullfile(figures_dir, 'Statistics_Complete_v2.csv'));
fprintf('  ✓ Saved Statistics_Complete_v2.csv\n\n');

fprintf('=============================================\n');
fprintf('ADDITIONAL ANALYSES COMPLETE\n');
fprintf('New figures: FigS1, FigS2, FigS3\n');
fprintf('Updated stats: Statistics_Complete_v2.csv\n');
fprintf('=============================================\n');%% ====================================================================
%  HELPER: Benjamini-Hochberg FDR correction
%  ====================================================================
function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q)
if nargin < 2, q = 0.05; end
pvals = pvals(:);
n     = length(pvals);
[sorted_p, sort_idx] = sort(pvals);
ecdf_factor = (1:n)' / n;
below  = sorted_p <= q * ecdf_factor;
crit_p = 0;
if any(below), crit_p = max(sorted_p(below)); end
h = pvals <= crit_p;

% Adjusted p-values
adj_p_sorted = nan(n,1);
for i = n:-1:1
    if i == n
        adj_p_sorted(i) = sorted_p(i);
    else
        adj_p_sorted(i) = min(sorted_p(i) * n / i, adj_p_sorted(i+1));
    end
end
adj_p_sorted = min(adj_p_sorted, 1);
adj_p = nan(n,1);
adj_p(sort_idx) = adj_p_sorted;
adj_ci_cvrg = 1 - q;
end
