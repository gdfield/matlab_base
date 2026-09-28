function os_make_figures_claude_auto(which_fig)
% os_make_figures_claude_auto  QC / results figures for the claude_auto analyses.
%   os_make_figures_claude_auto('pairing' | 'classifier' | 'axons' | 'orientation' | 'cluster20121031' | 'framerate')
% Figures are drawn in hidden windows and written as PDFs to os_analysis/claude_figures/.
% 2026-09-27 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
fd = [oa 'claude_figures/'];
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 8.5], 'DefaultAxesFontSize', 8, 'DefaultTextFontSize', 8, 'DefaultLegendFontSize', 7, 'DefaultAxesTitleFontWeight', 'normal');
switch which_fig
    case 'pairing',         fig_pairing(oa, f);
    case 'classifier',      fig_classifier(oa, f);
    case 'axons',           fig_axons(oa, f);
    case 'orientation',     fig_orientation(oa, f);
    case 'cluster20121031', fig_cluster(oa, f);
    case 'framerate',       fig_framerate(oa, f);
    case 'types_mi',        fig_types_mi(oa, f);
    case 'types_rf',        fig_types_rf(oa, f);
    case 'old_compare',     fig_old_compare(oa, f);
end
exportgraphics(f, [fd, 'fig_', which_fig, '.pdf'], 'ContentType', 'vector');
close(f);
end

function fig_pairing(oa, f)
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
S = load([oa 'pairing_shift_test_claude_auto.mat'], 'PS'); PS = S.PS; n = numel(PS);
[~, o] = sort({PS.dataset}); PS = PS(o);
ax = nexttile(tl, 1); hold(ax, 'on');
for i = 1:n
    plot(ax, [i i], [PS(i).null_mean PS(i).null_max], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 4);
    plot(ax, i, PS(i).max_other_shift, 's', 'Color', [0.85 0.4 0.1], 'MarkerFaceColor', [0.85 0.4 0.1]);
    plot(ax, i, PS(i).observed, 'o', 'Color', [0 0.3 0.7], 'MarkerFaceColor', [0 0.3 0.7]);
end
set(ax, 'XTick', 1:n, 'XTickLabel', {PS.dataset}, 'XTickLabelRotation', 60, 'TickDir', 'out'); xlim(ax, [0 n + 1]);
ylabel(ax, 'median within-condition reliability');
title(ax, 'Pairing check: observed (blue), best shifted (orange), shuffle null mean-max (grey)');
ax2 = nexttile(tl, 2); hold(ax2, 'on');
for i = 1:n, plot(ax2, PS(i).shifts, PS(i).shift_profile, '-', 'Color', [0 0.3 0.7 0.35]); end
xlabel(ax2, 'label shift k (trial t paired with file trial t+k)'); ylabel(ax2, 'median reliability');
title(ax2, 'Shift profiles, all 20 datasets: peak at k = 0 in every dataset'); set(ax2, 'TickDir', 'out');
end

function fig_classifier(oa, f)
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
S = load([oa 'claude_features/curation_classifier_results.mat'], 'O', 'N', 'R'); O = S.O; N = S.N; R = S.R;
ydef = O.label == 2; yos = O.label >= 1;
ax = nexttile(tl, 1); hold(ax, 'on');
[x1, y1] = perfcurve(ydef, O.P_def_lodo, true); plot(ax, x1, y1, 'LineWidth', 1.5);
[x2, y2] = perfcurve(yos, O.P_os_lodo, true); plot(ax, x2, y2, 'LineWidth', 1.5); plot(ax, [0 1], [0 1], 'k:');
legend(ax, sprintf('definite vs rest, AUC %.2f', R.auc_def_full), sprintf('OS vs rejected, AUC %.2f', R.auc_os_full), 'Location', 'southeast');
xlabel(ax, 'false positive rate'); ylabel(ax, 'true positive rate'); title(ax, 'Leave-one-dataset-out ROC'); axis(ax, 'square');
feats = {'F2n_mean', 'frac_cond_os', 'split_half_best', 'osi_best', 'corr_best'};
lab = categorical(O.label, [0 1 2], {'rejected', 'maybe', 'definite'});
for k = 1:4
    ax = nexttile(tl, k + 1);
    boxchart(ax, lab, O.(feats{k}), 'MarkerStyle', '.'); title(ax, strrep(feats{k}, '_', '\_')); set(ax, 'TickDir', 'out');
end
ax = nexttile(tl, 6); hold(ax, 'on');
e = 0:0.05:1;
histogram(ax, O.P_def_lodo(O.label == 0), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2);
histogram(ax, O.P_def_lodo(O.label == 1), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2);
histogram(ax, O.P_def_lodo(O.label == 2), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2);
histogram(ax, N.P_def, e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'k');
xline(ax, R.t_def, '--', 'definite cutoff'); legend(ax, 'old rejected', 'old maybe', 'old definite', 'new candidates', 'Location', 'north');
xlabel(ax, 'P(definite)'); ylabel(ax, 'fraction'); title(ax, 'P(definite): old (LODO) vs new');
end

function fig_axons(oa, f)
tl = tiledlayout(f, 6, 10, 'TileSpacing', 'tight', 'Padding', 'compact');
S = load([oa 'claude_axons/population_axon_axes.mat'], 'P'); P = S.P; n = height(P);
nr = 6; nc = ceil(n / nr);
for i = 1:n
    ax = polaraxes(tl); ax.Layout.Tile = i;
    A = load([oa 'claude_axons/' char(P.dataset(i)) '_axons.mat']);
    if isfield(A, 'AX'), a = A.AX.ax_angle(A.AX.good); polarhistogram(ax, deg2rad(a), 24, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); end
    ax.ThetaTickLabel = []; ax.RTickLabel = [];
    title(ax, {char(P.dataset(i)), sprintf('R=%.2f', P.dir_R(i))}, 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, 'EI axon propagation directions (array coordinates), cells with a clear axon; R = resultant length', 'FontSize', 9);
end

function fig_orientation(oa, f)
tl = tiledlayout(f, 6, 10, 'TileSpacing', 'tight', 'Padding', 'compact');
S = load([oa 'claude_features/feature_table.mat'], 'T'); T = S.T;
A = load([oa 'claude_axons/population_axon_axes.mat'], 'P'); P = A.P;
ds = unique(T.dataset(T.label == 2 | T.is_new)); n = numel(ds); nr = 6; nc = ceil(n / nr);
for i = 1:n
    ax = polaraxes(tl); ax.Layout.Tile = i; hold(ax, 'on');
    s = strcmp(T.dataset, ds{i}) & (T.label == 2 | T.is_new);
    psi = T.pref_ori_best(s); psi = [psi; psi + 180];
    polarhistogram(ax, deg2rad(psi), 36, 'FaceColor', [0 0.3 0.7], 'EdgeColor', 'none');
    j = find(P.dataset == ds{i});
    if ~isempty(j), r = ax.RLim(2); polarplot(ax, deg2rad(P.axis(j) + [0 180]), [r r], '-', 'Color', [0.85 0.4 0.1], 'LineWidth', 1.5); end
    ax.ThetaTickLabel = []; ax.RTickLabel = [];
    lb = 'old def'; if any(T.is_new(s)), lb = 'new cand'; end
    title(ax, {ds{i}, sprintf('%s n=%d', lb, sum(s))}, 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, 'Preferred drift axis (blue, stimulus coords) vs population axon axis (orange, array coords; frames not aligned)', 'FontSize', 9);
end

function fig_cluster(oa, f)
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
S = load([oa 'claude_features/d20121031_definite_clusters.mat'], 'D1'); D1 = S.D1;
isA = D1.cluster == 'A'; isH = D1.vtype == "OFF h os";
ax = polaraxes(tl); ax.Layout.Tile = 1; hold(ax, 'on');
polarhistogram(ax, deg2rad([D1.pref_ori_best(isA); D1.pref_ori_best(isA) + 180]), 36, 'FaceColor', [0.85 0.4 0.1], 'EdgeColor', 'none');
polarhistogram(ax, deg2rad([D1.pref_ori_best(~isA); D1.pref_ori_best(~isA) + 180]), 36, 'FaceColor', [0 0.3 0.7], 'EdgeColor', 'none');
title(ax, '2012-10-31-0 definite: drift axis (orange 81, blue 171 deg)');
ax = nexttile(tl, 2); hold(ax, 'on');
g = categorical(repmat("81-deg group", height(D1), 1)); g(~isA) = "171-deg group"; g(isH) = "171-deg, Vision hOS";
boxchart(ax, g, D1.mod_index); ylabel(ax, 'modulation index 2|F1|/F0'); title(ax, 'Simple/complex index by group'); set(ax, 'TickDir', 'out');
ax = nexttile(tl, 3); hold(ax, 'on');
scatter(ax, D1.pref_ori_best(isA), D1.rf_axis(isA), 25, [0.85 0.4 0.1], 'filled');
scatter(ax, D1.pref_ori_best(~isA), D1.rf_axis(~isA), 25, [0 0.3 0.7], 'filled');
xlabel(ax, 'preferred drift axis (deg, stimulus)'); ylabel(ax, 'RF long axis (deg, STA image)'); xlim(ax, [0 180]); ylim(ax, [0 180]); axis(ax, 'square');
title(ax, 'RF long axis vs drift axis'); set(ax, 'TickDir', 'out');
ax = nexttile(tl, 4); hold(ax, 'on');
boxchart(ax, categorical(cellstr(D1.cluster)), D1.rf_area); ylabel(ax, 'RF area (pixels, |z| > 4)'); title(ax, 'RF area (A = 81 deg, B = 171 deg)'); set(ax, 'TickDir', 'out');
end

function fig_framerate(oa, f)
tl = tiledlayout(f, 1, 1, 'Padding', 'compact');
d = dir([oa 'claude_modulation/*_mi.mat']); nm = {}; amp = [];
for k = 1:numel(d)
    S = load([d(k).folder '/' d(k).name]); if ~isfield(S, 'MI'), continue; end
    u = S.MI.Properties.UserData; nm{end + 1} = u.dataset; amp(end + 1, :) = mean(u.fr_amp, 2, 'omitnan')'; %#ok<AGROW>
end
[nm, o] = sort(nm); amp = amp(o, :);
ax = nexttile(tl); imagesc(ax, amp); colorbar(ax); colormap(ax, parula);
set(ax, 'XTick', 1:5, 'XTickLabel', {'30', '60', '75', '85', '120'}, 'YTick', 1:numel(nm), 'YTickLabel', nm, 'FontSize', 6);
xlabel(ax, 'candidate display frame rate (Hz)'); title(ax, 'Spectral SNR at TF = rate/TP (chosen = row max)');
end

function fig_types_mi(oa, f)
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); OS = S.OS;
tl = tiledlayout(f, 6, 10, 'TileSpacing', 'tight', 'Padding', 'compact');
ax = nexttile(tl, 1, [2 3]); hold(ax, 'on'); e = 0:0.1:2.2;
histogram(ax, OS.MI(OS.status == "old_definite"), e, 'DisplayStyle', 'stairs', 'LineWidth', 1.5);
histogram(ax, OS.MI(OS.status == "old_maybe"), e, 'DisplayStyle', 'stairs', 'LineWidth', 1.5);
histogram(ax, OS.MI(OS.is_new), e, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'k');
xline(ax, 0.7, '--'); legend(ax, 'old definite', 'old maybe', 'new def+maybe', 'Location', 'northwest');
xlabel(ax, 'modulation index 2|F1|/F0 at TF'); ylabel(ax, 'cells'); title(ax, 'Simple (right) vs complex (left) OS cells');
ds = unique(OS.dataset); k = 3;
for i = 1:numel(ds)
    s = strcmp(OS.dataset, ds{i}) & ~isnan(OS.MI); if sum(s) < 5, continue; end
    k = k + 1; if mod(k - 1, 10) < 3 && k <= 23, k = k + (3 - mod(k - 1, 10)); end
    if k > 60, break; end
    ax = polaraxes(tl); ax.Layout.Tile = k; hold(ax, 'on');
    a = OS.pref_axis(s); cx = OS.cx(s);
    r = 0.3 + 0.7 * rand(sum(s), 1);
    polarscatter(ax, deg2rad([a(~cx); a(~cx) + 180]), [r(~cx); r(~cx)], 6, [0 0.3 0.7], 'filled');
    polarscatter(ax, deg2rad([a(cx); a(cx) + 180]), [r(cx); r(cx)], 10, [0.85 0.1 0.1], 'filled');
    ax.ThetaTickLabel = []; ax.RTickLabel = []; rlim(ax, [0 1]);
    title(ax, {ds{i}, sprintf('S%d C%d', sum(~cx), sum(cx))}, 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, 'Preferred drift axis per dataset: simple (blue, MI>=0.7) vs complex (red, MI<0.7); definite + maybe', 'FontSize', 9);
end

function fig_types_rf(oa, f)
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); OS = S.OS; OS = OS(~isnan(OS.rf_area), :);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
g = categorical(repmat("simple", height(OS), 1)); g(OS.cx) = "complex";
v = {'rf_area_rel', 'rf_aspect', 'rf_frac_largest', 'rf_npatch', 'rf_peakz'};
lb = {'RF area / dataset median', 'RF elongation (aspect)', 'fraction of RF in largest patch', 'number of RF patches', 'STA peak |z|'};
for k = 1:5
    ax = nexttile(tl, k); boxchart(ax, g, OS.(v{k}), 'MarkerStyle', '.'); ylabel(ax, lb{k}); set(ax, 'TickDir', 'out');
    p = ranksum(OS.(v{k})(OS.cx), OS.(v{k})(~OS.cx)); title(ax, sprintf('ranksum p = %.2g (n S=%d, C=%d)', p, sum(~OS.cx), sum(OS.cx)));
end
ax = nexttile(tl, 6); hold(ax, 'on');
histogram(ax, OS.rf_minus_pref(~OS.cx & OS.rf_aspect > 1.5), -90:10:90, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.5);
histogram(ax, OS.rf_minus_pref(OS.cx & OS.rf_aspect > 1.5), -90:10:90, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.5);
legend(ax, 'simple', 'complex'); xlabel(ax, 'RF long axis - preferred drift axis (deg)'); title(ax, 'Does RF elongation predict orientation? (aspect > 1.5)');
end

function fig_old_compare(oa, f)
S = load([oa 'claude_features/old_datasets_claude_vs_gdf_classification.mat'], 'C'); C = S.C;
P = readtable([oa 'claude_features/old_datasets_comparison_per_dataset.csv']);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, 1); [tab, ~, ~, lb] = crosstab(C.gdf_label, C.claude_label);
imagesc(ax, tab ./ sum(tab, 2)); colormap(ax, flipud(gray)); clim(ax, [0 1]); colorbar(ax);
for a = 1:size(tab, 1), for b = 1:size(tab, 2), text(ax, b, a, sprintf('%d', tab(a, b)), 'HorizontalAlignment', 'center', 'Color', [0.85 0.3 0.1], 'FontWeight', 'bold'); end, end
set(ax, 'XTick', 1:4, 'XTickLabel', strrep(categories(C.claude_label), '_', ' '), 'YTick', 1:3, 'YTickLabel', categories(C.gdf_label));
xlabel(ax, 'Claude (leave-one-dataset-out)'); ylabel(ax, 'GDF curation'); title(ax, 'Counts; shading = fraction of each GDF row');
ax = nexttile(tl, 2, [1 2]); P = sortrows(P, 'dataset');
bar(ax, [P.gdf_def P.claude_def], 'grouped'); set(ax, 'XTick', 1:height(P), 'XTickLabel', P.dataset, 'XTickLabelRotation', 70, 'FontSize', 6, 'TickDir', 'out');
legend(ax, 'GDF definite', 'Claude definite'); ylabel(ax, 'cells'); title(ax, 'Definite OS cells per dataset');
ax = nexttile(tl, 5, [1 2]);
bar(ax, [P.gdf_def + P.gdf_maybe, P.claude_def + P.claude_maybe], 'grouped'); set(ax, 'XTick', 1:height(P), 'XTickLabel', P.dataset, 'XTickLabelRotation', 70, 'FontSize', 6, 'TickDir', 'out');
legend(ax, 'GDF definite + maybe', 'Claude definite + maybe'); ylabel(ax, 'cells'); title(ax, 'Definite + maybe per dataset (GDF count includes cells failing the current gate)');
ax = nexttile(tl, 4); c = C.claude_label ~= "not_candidate"; hold(ax, 'on');
e = 0:0.05:1; cols = lines(3); gl = categories(C.gdf_label);
for k = 1:3, histogram(ax, C.P_def_lodo(c & C.gdf_label == gl{k}), e, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', cols(k, :), 'Normalization', 'probability'); end
xline(ax, 0.73, '--'); legend(ax, strcat('GDF ', gl), 'Location', 'north'); xlabel(ax, 'P(definite), leave-one-dataset-out'); ylabel(ax, 'fraction');
title(ax, 'Model confidence by your label');
end
