function os_robust_explainer_claude_auto(dataset, cell_ids)
% os_robust_explainer_claude_auto  Figures explaining the bias-corrected bootstrap
% lower bound (r2_lbx, M10) with example cells, plus population summaries.
% fig_robust_explainer.pdf: one row per example cell, best condition (from
%   robust_cells.mat): (1) direction tuning, mean count per direction with the
%   5-95% bootstrap range; (2) distributions of r2 over bootstrap replicates (blue)
%   and over the permutation null (grey), with observed r2, its 5th percentile
%   (r2_lb), the null mean (noise floor) and r2_lbx = r2_lb - null mean.
% fig_robust_population.pdf: noise floor vs spike count; r2_lbx by label;
%   r2_lbx vs reliability; per-dataset fraction of GDF definite cells above 0.16.
% Recomputes the example cells with os_robust_tuning_claude_auto (same seed, so
% values equal the batch values; checked and printed).
% 2026-09-28 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
fd = [oa 'claude_figures/robust/'];
L = load([oa 'claude_robust/robust_cells.mat'], 'P'); P = L.P;
reg = load([oa 'claude_features/registry_used.mat']); dl = reg.allds(strcmp(reg.nm, dataset));
evalc('R = os_robust_tuning_claude_auto(dl, ''save'', false, ''keep_ids'', cell_ids);');
n = numel(cell_ids);
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 10 2.6 * n], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, n, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:n
    pr = P(P.dataset == dataset & P.cell_id == cell_ids(i), :);
    j = find(R.cond_sp == pr.best_sp & R.cond_tp == pr.best_tp, 1);
    kc = find(R.D.kc_ids == cell_ids(i)); r = find(R.cell_ids == cell_ids(i));
    Cc = R.D.Cc{j}(kc, :); dk = R.D.dk{j}; cnt = R.D.boot_cnt{j}; nrep = R.D.nrep{j}; nd = numel(R.dirs);
    m = accumarray(dk(:), Cc(:), [nd 1]) ./ nrep(:);
    mb = zeros(nd, size(cnt, 2));
    for k = 1:nd, t = dk == k; mb(k, :) = (Cc(t) * cnt(t, :)) / nrep(k); end
    lo = prctile(mb, 5, 2); hi = prctile(mb, 95, 2);
    br = R.D.boot_r2{j}(kc, :); nr = R.D.null_r2{j}(kc, :);
    r2 = R.r2(r, j); lb = R.r2_lb(r, j); nm = R.r2_null(r, j);
    fprintf('cell %d: r2 %.4f (table %.4f), lb %.4f (table %.4f), lbx %.4f (table %.4f)\n', cell_ids(i), r2, pr.r2_best, lb, pr.r2_lb_best, lb - nm, pr.r2_lbx_best);
    ax = nexttile(tl); hold(ax, 'on'); x = R.dirs(:);
    fill(ax, [x; flipud(x)], [lo; flipud(hi)], [0.75 0.85 1], 'EdgeColor', 'none');
    plot(ax, x, m, 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 3); xlim(ax, [min(x) - 10, max(x) + 10]);
    xticks(ax, x); xlabel(ax, 'drift direction (deg)'); ylabel(ax, 'mean spikes / trial');
    title(ax, sprintf('%s cell %d (%s)\\newlineSP %g TP %g, %d spikes, %d reps', dataset, cell_ids(i), pr.label, pr.best_sp, pr.best_tp, pr.nspk_best, min(nrep)), 'FontSize', 7);
    ax = nexttile(tl); hold(ax, 'on'); e = linspace(0, max([br(:); nr(:); r2]) * 1.05 + eps, 60);
    histogram(ax, nr, e, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'probability');
    histogram(ax, br, e, 'FaceColor', [0.2 0.45 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'Normalization', 'probability');
    yl = ylim(ax); xline(ax, r2, 'k-', 'LineWidth', 1.3); xline(ax, lb, 'b--', 'LineWidth', 1.3); xline(ax, nm, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.3);
    plot(ax, [nm lb], 0.9 * yl(2) * [1 1], 'r-', 'LineWidth', 2);
    xlabel(ax, 'r2 (vector-sum orientation index)'); ylabel(ax, 'fraction');
    legend(ax, {'permutation null', 'bootstrap', 'observed r2', '5th pct (r2 lb)', 'null mean', 'r2 lbx'}, 'Box', 'off', 'FontSize', 6, 'Location', 'best');
    title(ax, sprintf('r2 %.2f, lb %.2f, null mean %.2f  ->  r2 lbx %.2f', r2, lb, nm, lb - nm), 'FontSize', 7);
    ax = nexttile(tl); axis(ax, 'off');
    text(ax, 0, 0.9, sprintf(['q (FDR) %.2g\nreliability across reps %.2f\nr1 (direction) %.2f\n' ...
        'conditions significant %d of %d\npreferred axis %.0f deg, bootstrap SD %.1f deg'], ...
        pr.q_min, pr.rel_best, pr.r1_best, pr.k_sig, pr.n_cond, pr.pref_axis_best, pr.axis_sd_best), 'FontSize', 7, 'VerticalAlignment', 'top');
end
title(tl, 'How the bias-corrected bootstrap lower bound (r2 lbx, red bar) is formed', 'FontSize', 9);
exportgraphics(f, [fd 'fig_robust_explainer.pdf'], 'ContentType', 'vector'); close(f);
% ---------- population
old = ismember(P.label, ["definite", "maybe", "gate_rejected", "background"]);
lab = ["definite", "maybe", "gate_rejected", "background"]; cl = [0.1 0.35 0.8; 0.3 0.7 0.9; 0.85 0.4 0.1; 0.6 0.6 0.6];
rng(3); sub = rand(height(P), 1) < 0.1;
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 10 7], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
for i = [4 3 2 1], m = old & P.label == lab(i) & (i < 4 | sub); scatter(ax, P.nspk_best(m), P.r2_best(m) - P.r2x_best(m), 5, cl(i, :), 'filled', 'MarkerFaceAlpha', 0.4); end
set(ax, 'XScale', 'log'); xlabel(ax, 'spikes in best condition'); ylabel(ax, 'null mean r2 (noise floor)');
title(ax, 'the r2 expected from noise alone falls with spike count'); legend(ax, {'background (10%)', 'gate rejected', 'maybe', 'definite'}, 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); e = -0.1:0.02:0.8;
for i = 1:4, m = old & P.label == lab(i); histogram(ax, P.r2_lbx_best(m), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', cl(i, :), 'LineWidth', 1.3); end
xline(ax, 0.16, 'k:'); xlabel(ax, 'r2 lbx (best condition)'); ylabel(ax, 'fraction of label');
legend(ax, cellstr(lab), 'Box', 'off', 'Interpreter', 'none'); title(ax, 'bias-corrected lower bound by label (earlier datasets)');
ax = nexttile(tl); hold(ax, 'on');
for i = [4 3 2 1], m = old & P.label == lab(i) & (i < 4 | sub); scatter(ax, P.rel_best(m), P.r2_lbx_best(m), 5, cl(i, :), 'filled', 'MarkerFaceAlpha', 0.4); end
yline(ax, 0.16, 'k:'); xlabel(ax, 'reliability across repetitions (best condition)'); ylabel(ax, 'r2 lbx'); title(ax, 'r2 lbx vs reliability');
ax = nexttile(tl); ds = unique(P.dataset(P.label == "definite")); fr = nan(numel(ds), 1); nn = fr;
for k = 1:numel(ds), m = P.dataset == ds(k) & P.label == "definite"; fr(k) = mean(P.r2_lbx_best(m) >= 0.16); nn(k) = sum(m); end
bar(ax, fr); set(ax, 'XTick', 1:numel(ds), 'XTickLabel', ds, 'XTickLabelRotation', 70, 'TickLabelInterpreter', 'none');
ylabel(ax, 'fraction of GDF definite with r2 lbx >= 0.16'); title(ax, 'per dataset (n definite shown above bars)');
text(ax, 1:numel(ds), fr + 0.03, string(nn), 'HorizontalAlignment', 'center', 'FontSize', 5); ylim(ax, [0 1.1]);
exportgraphics(f, [fd 'fig_robust_population.pdf'], 'ContentType', 'vector'); close(f);
end
