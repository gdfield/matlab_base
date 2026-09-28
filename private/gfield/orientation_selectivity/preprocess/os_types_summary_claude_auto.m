function Sm = os_types_summary_claude_auto()
% os_types_summary_claude_auto  Cross-dataset summary of os_types_hsimple_claude_auto:
% per-dataset anchor axis, class counts, ambiguity flag, RF availability, and RF
% polarity composition (STA time course, guess_polarity) of h-simple and complex
% vOS cells. Writes claude_types/summary_all_datasets.csv and
% claude_types/fig_types_summary_all.pdf. 2026-09-27 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/claude_types/';
d = dir([oa '*/types.mat']); n = numel(d); z = nan(n, 1);
Sm = table(strings(n, 1), z, z, z, z, z, false(n, 1), false(n, 1), z, z, z, z, z, strings(n, 1), z, z, false(n, 1), z, 'VariableNames', ...
    {'dataset', 'anchor_axis', 'n_h_simple', 'n_complex_vOS', 'n_orth_simple', 'n_other', 'ambiguous', 'has_rf', ...
    'h_on', 'h_off', 'cx_on', 'cx_off', 'cx_by_MI', 'anchor_method', 'anchor_margin', 'anchor_p_sep', 'excluded', 'h_med_flank'});
for k = 1:n
    L = load(fullfile(d(k).folder, d(k).name), 'T', 'info'); T = L.T; I = L.info;
    Sm.dataset(k) = I.dataset; Sm.anchor_axis(k) = I.anchor_axis; Sm.n_h_simple(k) = I.n_h_simple;
    Sm.n_complex_vOS(k) = I.n_complex_vOS; Sm.n_orth_simple(k) = I.n_orth_simple; Sm.n_other(k) = I.n_other;
    Sm.ambiguous(k) = I.ambiguous_anchor; Sm.has_rf(k) = I.has_rf;
    h = T.class == "h_simple"; c = T.class == "complex_vOS";
    Sm.h_on(k) = sum(h & T.polarity > 0); Sm.h_off(k) = sum(h & T.polarity < 0);
    Sm.cx_on(k) = sum(c & T.polarity > 0); Sm.cx_off(k) = sum(c & T.polarity < 0);
    Sm.cx_by_MI(k) = sum(c & contains(T.reason, "MI"));
    Sm.anchor_method(k) = I.anchor_method; Sm.anchor_margin(k) = I.anchor_margin; Sm.anchor_p_sep(k) = I.anchor_p_sep;
    Sm.excluded(k) = I.excluded; Sm.h_med_flank(k) = median(T.flank_score(h), 'omitnan');
end
Sm = sortrows(Sm, 'dataset'); writetable(Sm, [oa 'summary_all_datasets.csv']);
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 7], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); b = bar(ax, [Sm.n_h_simple Sm.n_complex_vOS Sm.n_orth_simple Sm.n_other], 'stacked');
cl = [0.1 0.35 0.8; 0.85 0.2 0.15; 0.2 0.6 0.3; 0.7 0.7 0.7]; for i = 1:4, b(i).FaceColor = cl(i, :); end
set(ax, 'XTick', 1:n, 'XTickLabel', Sm.dataset, 'XTickLabelRotation', 70, 'TickLabelInterpreter', 'none'); ylabel(ax, 'cells');
hold(ax, 'on'); a = find(Sm.excluded); plot(ax, a, sum(table2array(Sm(a, 3:6)), 2) + 3, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
nr = find(~Sm.has_rf); plot(ax, nr, sum(table2array(Sm(nr, 3:6)), 2) + 7, 'ko', 'MarkerSize', 4);
legend(ax, {'h-simple', 'complex vOS', 'orth-simple', 'other', 'excluded (<5 h-simple)', 'no WN RFs (count-based anchor)'}, 'Box', 'off', 'Location', 'northwest');
title(ax, 'Class counts per dataset (curated OS cells, def + maybe)');
ax = nexttile(tl); hold(ax, 'on'); w = 0.38; x = 1:n;
bar(ax, x - w / 2, [Sm.h_off Sm.h_on], w, 'stacked'); bb = bar(ax, x + w / 2, [Sm.cx_off Sm.cx_on], w, 'stacked');
ch = ax.Children; set(ch(end), 'FaceColor', [0.05 0.2 0.55]); set(ch(end - 1), 'FaceColor', [0.55 0.7 0.95]);
set(bb(1), 'FaceColor', [0.6 0.1 0.08]); set(bb(2), 'FaceColor', [0.98 0.65 0.6]);
set(ax, 'XTick', 1:n, 'XTickLabel', Sm.dataset, 'XTickLabelRotation', 70, 'TickLabelInterpreter', 'none'); ylabel(ax, 'cells with RF polarity');
legend(ax, {'h-simple OFF', 'h-simple ON', 'complex vOS OFF', 'complex vOS ON'}, 'Box', 'off', 'Location', 'northwest');
title(ax, sprintf('RF polarity (STA time course): h-simple %d OFF / %d ON; complex vOS %d OFF / %d ON', ...
    sum(Sm.h_off), sum(Sm.h_on), sum(Sm.cx_off), sum(Sm.cx_on)));
exportgraphics(f, [oa 'fig_types_summary_all.pdf'], 'ContentType', 'vector'); close(f);
end
