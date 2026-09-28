function A = os_robust_analysis_claude_auto()
% os_robust_analysis_claude_auto  Figures and statistics for the trial-level
% classification metrics (os_robust_tuning_claude_auto / os_robust_compare_claude_auto).
% Earlier (GDF-curated) datasets only unless stated. Writes
% os_analysis/claude_figures/robust/fig_robust_*.pdf and claude_robust/robust_auc.csv.
%   calibration : p-value histograms, real vs trials shuffled across directions
%                 (2012-10-31-0, 2022-12-21-0)
%   bootstrap   : r2 vs bootstrap lower bound, colored by spike count
%   roc         : ROC, GDF definite vs gate-rejected, per metric; AUC by 8/12 dirs
%   sampling    : r2 and old OSI of GDF definite cells, 8- vs 12-direction datasets
%   condspec    : fraction of conditions significant, by label; M8 model-miss cells
%   rule        : sweep of rule (k_sig >= 1 & r2_lb_best >= x): fraction flagged per label
% 2026-09-28 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
fd = [oa 'claude_figures/robust/']; if ~exist(fd, 'dir'), mkdir(fd); end
L = load([oa 'claude_robust/robust_cells.mat'], 'P'); P = L.P;
old = ismember(P.label, ["definite", "maybe", "gate_rejected", "background"]);
lab = ["definite", "maybe", "gate_rejected", "background"]; cl = [0.1 0.35 0.8; 0.3 0.7 0.9; 0.85 0.4 0.1; 0.6 0.6 0.6];
A = struct();
% ---------- calibration
reg = load([oa 'claude_features/registry_used.mat']);
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.2], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); dsc = {'2012-10-31-0', '2022-12-21-0'};
for k = 1:2
    dl = reg.allds(strcmp(reg.nm, dsc{k}));
    evalc('Rr = os_robust_tuning_claude_auto(dl, ''save'', false);');
    evalc('Rc = os_robust_tuning_claude_auto(dl, ''calib_shuffle'', true, ''save'', false, ''seed'', 7);');
    q = Rr.nspk >= 20; ax = nexttile(tl); hold(ax, 'on'); e = 0:0.025:1;
    histogram(ax, Rr.p_r2(q), e, 'Normalization', 'probability', 'FaceColor', [0.2 0.4 0.8]);
    histogram(ax, Rc.p_r2(q), e, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
    A.calib.(['d' strrep(dsc{k}, '-', '_')]) = [mean(Rc.p_r2(q) < 0.05), mean(Rr.p_r2(q) < 0.05)];
    legend(ax, {sprintf('real (p<0.05: %.2f)', mean(Rr.p_r2(q) < 0.05)), sprintf('directions shuffled (p<0.05: %.3f)', mean(Rc.p_r2(q) < 0.05))}, 'Box', 'off');
    xlabel(ax, 'permutation p (r2)'); ylabel(ax, 'fraction of cell x condition tests'); title(ax, sprintf('%s, all cells with >= 20 spikes', dsc{k}));
end
exportgraphics(f, [fd 'fig_robust_calibration.pdf'], 'ContentType', 'vector'); close(f);
% ---------- bootstrap
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.6], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); s = old & P.gate_pass | ismember(P.label, ["definite", "maybe"]);
scatter(ax, P.r2_best(s), P.r2_lbx_best(s), 8, log10(P.nspk_best(s) + 1), 'filled', 'MarkerFaceAlpha', 0.6); colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = 'log10 spikes (best condition)'; hold(ax, 'on'); plot(ax, [0 1], [0 1], 'k:');
xlabel(ax, 'r2 (vector-sum orientation index)'); ylabel(ax, 'bootstrap 5th pct - null mean (r2 lbx)'); title(ax, 'curated + gate-passing cells, earlier datasets');
ax = nexttile(tl); hold(ax, 'on');
for i = 1:4, m = old & P.label == lab(i); if i == 4, m = m & rand(height(P), 1) < 0.1; end
    scatter(ax, log10(P.nspk_best(m) + 1), P.r2_best(m) - P.r2x_best(m), 6, cl(i, :), 'filled', 'MarkerFaceAlpha', 0.4); end
xlabel(ax, 'log10 spikes (best condition)'); ylabel(ax, 'null mean r2 (noise floor)'); legend(ax, [cellstr(lab(1:3)) {'background (10% sample)'}], 'Box', 'off', 'Interpreter', 'none');
title(ax, 'noise floor of r2 rises as spike count falls (corrected for)');
exportgraphics(f, [fd 'fig_robust_bootstrap.pdf'], 'ContentType', 'vector'); close(f);
% ---------- roc + AUC table
mets = {'r2_lbx_best', 'r2x_best', 'r2_lb_best', 'q_min', 'frac_sig', 'rel_best', 'osi_best_old', 'P_def_lodo'};
nm = {'r2 lower bound, bias-corr.', 'r2 bias-corr.', 'r2 lower bound (raw)', '-log10 q', 'frac. conditions sig.', 'reliability', 'old OSI (gate)', 'M5 P(def) (trained on these labels)'};
cmp = {{"definite_gp", "gate_rejected"}, {["definite_gp", "maybe_gp"], "gate_rejected"}, {"definite", "background"}};
P.label2 = P.label; P.label2(P.label == "definite" & P.gate_pass) = "definite_gp"; P.label2(P.label == "maybe" & P.gate_pass) = "maybe_gp";
P.label2(P.label == "definite" & ~P.gate_pass) = "definite_ngp"; P.label2(P.label == "maybe" & ~P.gate_pass) = "maybe_ngp";
rows = {};
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 3.6], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); cm = lines(numel(mets));
for c = 1:3
    pos = ismember(P.label2, cmp{c}{1}) | (c == 3 & P.label == "definite"); neg = ismember(P.label2, cmp{c}{2});
    ax = nexttile(tl); hold(ax, 'on'); lg = {};
    for m = 1:numel(mets)
        x = P.(mets{m}); if strcmp(mets{m}, 'q_min'), x = -log10(x); end
        for nd_ = [0 8 12]
            s = (pos | neg) & ~isnan(x); if nd_ > 0, s = s & P.n_dirs == nd_; end
            if sum(s & pos) < 10 || sum(s & neg) < 10, a = NaN; X = []; Y = [];
            else, [X, Y, ~, a] = perfcurve(pos(s), x(s), true); end
            rows(end + 1, :) = {strjoin(cmp{c}{1}, '+') + " vs " + cmp{c}{2}, mets{m}, nd_, sum(s & pos), sum(s & neg), a}; %#ok<AGROW>
            if nd_ == 0 && ~isempty(X), plot(ax, X, Y, 'Color', cm(m, :), 'LineWidth', 1); lg{end + 1} = sprintf('%s %.2f', nm{m}, a); end %#ok<AGROW>
        end
    end
    plot(ax, [0 1], [0 1], 'k:'); legend(ax, lg, 'Box', 'off', 'Location', 'southeast', 'FontSize', 6);
    xlabel(ax, 'false positive rate'); ylabel(ax, 'true positive rate');
    title(ax, sprintf('%s vs %s (n %d / %d)', strjoin(cmp{c}{1}, '+'), cmp{c}{2}, sum(pos), sum(neg)), 'Interpreter', 'none');
end
exportgraphics(f, [fd 'fig_robust_roc.pdf'], 'ContentType', 'vector'); close(f);
Ta = cell2table(rows, 'VariableNames', {'comparison', 'metric', 'n_dirs_subset', 'n_pos', 'n_neg', 'AUC'});
writetable(Ta, [oa 'claude_robust/robust_auc.csv']); A.auc = Ta;
% ---------- sampling
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.2], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
vv = {'r2x_best', 'r2_lbx_best', 'osi_best_old'}; ttl = {'r2 bias-corr.', 'r2 lower bound bias-corr.', 'old OSI'};
for i = 1:3
    ax = nexttile(tl); hold(ax, 'on'); e = linspace(0, 1, 26);
    for nd_ = [8 12]
        m = P.label == "definite" & P.n_dirs == nd_; x = P.(vv{i})(m);
        histogram(ax, x, e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2);
        A.sampling.(vv{i})(nd_ / 4 - 1) = median(x, 'omitnan');
    end
    x8 = P.(vv{i})(P.label == "definite" & P.n_dirs == 8); x12 = P.(vv{i})(P.label == "definite" & P.n_dirs == 12);
    pr = ranksum(x8, x12);
    legend(ax, {sprintf('8 dirs, median %.2f', median(x8, 'omitnan')), sprintf('12 dirs, median %.2f', median(x12, 'omitnan'))}, 'Box', 'off');
    title(ax, sprintf('%s: 8 vs 12 dirs, p %.2g', ttl{i}, pr), 'FontSize', 7); xlabel(ax, [ttl{i} ' (GDF definite)']);
end
exportgraphics(f, [fd 'fig_robust_sampling.pdf'], 'ContentType', 'vector'); close(f);
% ---------- condition specificity
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.2], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on'); e = -0.05:0.1:1.05;
for i = 1:3, m = old & P.label == lab(i) & P.n_cond >= 3; histogram(ax, P.frac_sig(m), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2, 'EdgeColor', cl(i, :)); end
legend(ax, cellstr(lab(1:3)), 'Box', 'off', 'Interpreter', 'none'); xlabel(ax, 'fraction of conditions with significant orientation tuning (q<0.05, r2>r1)');
ylabel(ax, 'fraction of cells'); title(ax, 'datasets with >= 3 conditions');
ax = nexttile(tl); DL = load([oa 'claude_features/old_disagreements.mat']); fn = fieldnames(DL); DL = DL.(fn{1});
hold(ax, 'on'); m = P.label == "definite"; scatter(ax, P.frac_sig(m) + 0.02 * randn(sum(m), 1), P.r2_lbx_best(m), 8, [0.7 0.75 0.9], 'filled');
try
    km = string(DL.dataset) + "_" + string(DL.cell_id); sel = contains(string(DL.category), "not_OS") | contains(string(DL.category), "not");
    [tf, lo] = ismember(km(sel), P.dataset + "_" + string(P.cell_id)); lo = lo(tf);
    scatter(ax, P.frac_sig(lo) + 0.02 * randn(numel(lo), 1), P.r2_lbx_best(lo), 22, 'r', 'filled'); A.miss_cells = P(lo, {'dataset', 'cell_id', 'n_cond', 'k_sig', 'r2_best', 'r2_lbx_best', 'nspk_best', 'q_min'});
catch ME, A.miss_err = ME.message; end
xlabel(ax, 'fraction of conditions significant (jittered)'); ylabel(ax, 'r2 lower bound, bias-corrected (best condition)');
title(ax, sprintf('GDF definite (blue); red = %d M8 cells the M5 model called not OS', numel(lo)));
exportgraphics(f, [fd 'fig_robust_condspec.pdf'], 'ContentType', 'vector'); close(f);
% ---------- rule sweep
xs = 0:0.02:0.5; fr = zeros(numel(xs), 4); nf = fr;
for i = 1:numel(xs), for j = 1:4, m = old & P.label == lab(j); hit = P.k_sig(m) >= 1 & P.r2_lbx_best(m) >= xs(i); fr(i, j) = mean(hit); nf(i, j) = sum(hit); end, end
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.6], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on'); for j = 1:4, plot(ax, xs, fr(:, j), '-', 'Color', cl(j, :), 'LineWidth', 1.3); end
ylim(ax, [0 1]); xlabel(ax, 'threshold x on bias-corrected r2 lower bound'); ylabel(ax, 'fraction of label flagged');
legend(ax, cellstr(lab), 'Box', 'off', 'Interpreter', 'none', 'Location', 'northeast'); title(ax, 'rule: >= 1 significant condition (>= 50 spikes) AND r2 lbx >= x');
ax = nexttile(tl); hold(ax, 'on'); for j = 1:4, plot(ax, xs, nf(:, j), '-', 'Color', cl(j, :), 'LineWidth', 1.3); end
set(ax, 'YScale', 'log'); xlabel(ax, 'threshold x'); ylabel(ax, 'number of cells flagged'); title(ax, 'same, as counts (base rates differ: background n = 25,531)');
exportgraphics(f, [fd 'fig_robust_rule.pdf'], 'ContentType', 'vector'); close(f);
A.rule = array2table([xs(:) fr nf], 'VariableNames', ['x', cellstr(lab), strcat('n_', cellstr(lab))]);
A.gate_fail = groupsummary(P(ismember(P.label, ["definite", "maybe"]), :), {'label', 'gate_pass'}, 'median', {'r2_lbx_best', 'k_sig', 'rel_best', 'nspk_best'});
save([oa 'claude_robust/robust_analysis.mat'], 'A');
end
