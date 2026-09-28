function A = os_types_axis_plot_claude_auto()
% os_types_axis_plot_claude_auto  OSI vs MI vs preferred axis across datasets, with
% each dataset's preferred axes expressed relative to its h-simple anchor (M9 v3:
% 0 deg = mode of the h-simple group). Curated OS cells (def + maybe) of datasets
% kept by M9 (>= 5 h-simple). OSI = osi_max (pipeline OSI, max over conditions);
% MI = 2F1/F0. Note: h-simple cells sit near 0 by construction (they define it).
% A mirror in the newer rig could flip handedness between rigs, so the sign of the
% relative axis may not be comparable between 2012 and later datasets; the last
% panels therefore also show |relative axis| (0-90).
% Writes claude_figures/fig_types_axis_osi_mi.pdf and claude_types/axis_osi_mi_cells.csv.
% 2026-09-28 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
d = dir([oa 'claude_types/*/types.mat']); parts = {};
for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name), 'T', 'info'); I = L.info;
    if I.excluded || isnan(I.anchor_axis), continue; end
    T = L.T; t = table(T.dataset, T.cell_id, T.class, T.axis_rel_anchor, T.MI, T.osi_max, T.F2n, T.polarity, ...
        repmat(string(I.anchor_method), height(T), 1), 'VariableNames', {'dataset', 'cell_id', 'class', 'axis_rel', 'MI', 'OSI', 'F2n', 'polarity', 'anchor_method'});
    parts{end + 1} = t; %#ok<AGROW>
end
C = vertcat(parts{:}); C.era = repmat("later", height(C), 1); C.era(startsWith(string(C.dataset), "2012")) = "2012";
writetable(C, [oa 'claude_types/axis_osi_mi_cells.csv']);
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 7.5], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ok = ~isnan(C.axis_rel) & ~isnan(C.MI) & ~isnan(C.OSI);
ax = nexttile(tl); scatter(ax, C.axis_rel(ok), C.MI(ok), 10, C.OSI(ok), 'filled', 'MarkerFaceAlpha', 0.7); colormap(ax, turbo);
cb = colorbar(ax); cb.Label.String = 'OSI'; yline(ax, 0.7, 'k:'); xline(ax, [-22.5 22.5], 'k:'); xline(ax, [-67.5 67.5], 'k:');
xlim(ax, [-90 90]); xticks(ax, -90:45:90); xlabel(ax, 'preferred axis re h-simple anchor (deg)'); ylabel(ax, 'MI (F1/F0)');
title(ax, sprintf('all kept datasets (n = %d cells, %d datasets)', sum(ok), numel(unique(C.dataset))));
ax = nexttile(tl); hold(ax, 'on'); s = ok & C.MI >= 0.7; c = ok & C.MI < 0.7;
scatter(ax, C.axis_rel(s), C.OSI(s), 8, [0.1 0.35 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(ax, C.axis_rel(c), C.OSI(c), 12, [0.85 0.2 0.15], 'filled', 'MarkerFaceAlpha', 0.7);
xlim(ax, [-90 90]); xticks(ax, -90:45:90); xlabel(ax, 'preferred axis re anchor (deg)'); ylabel(ax, 'OSI');
legend(ax, {'MI >= 0.7 (simple)', 'MI < 0.7 (complex)'}, 'Box', 'off', 'Location', 'best'); title(ax, 'OSI vs axis');
ax = nexttile(tl); hold(ax, 'on'); e = -90:15:90;
histogram(ax, C.axis_rel(s), e, 'FaceColor', [0.1 0.35 0.8], 'FaceAlpha', 0.5);
histogram(ax, C.axis_rel(c), e, 'FaceColor', [0.85 0.2 0.15], 'FaceAlpha', 0.7);
xlim(ax, [-90 90]); xticks(ax, -90:45:90); xlabel(ax, 'preferred axis re anchor (deg)'); ylabel(ax, 'cells'); title(ax, 'axis distribution by MI class');
ax = nexttile(tl); scatter3(ax, C.axis_rel(ok), C.MI(ok), C.OSI(ok), 8, 1 + (C.MI(ok) < 0.7), 'filled'); colormap(ax, [0.1 0.35 0.8; 0.85 0.2 0.15]);
xlabel(ax, 'axis re anchor (deg)'); ylabel(ax, 'MI'); zlabel(ax, 'OSI'); xlim(ax, [-90 90]); view(ax, -35, 25); grid(ax, 'on'); title(ax, '3D: axis, MI, OSI');
ax = nexttile(tl); hold(ax, 'on'); mk = {'o', '^'}; er = ["2012", "later"];
for q = 1:2, m = ok & C.era == er(q); scatter(ax, abs(C.axis_rel(m)), C.MI(m), 12, C.OSI(m), 'filled', mk{q}, 'MarkerFaceAlpha', 0.7); end
colormap(ax, turbo); cb = colorbar(ax); cb.Label.String = 'OSI'; yline(ax, 0.7, 'k:'); xline(ax, [22.5 67.5], 'k:');
xlim(ax, [0 90]); xticks(ax, 0:15:90); xlabel(ax, '|preferred axis re anchor| (deg)'); ylabel(ax, 'MI');
legend(ax, {'2012 rig', 'later rig'}, 'Box', 'off', 'Location', 'best'); title(ax, 'folded (immune to a handedness flip)');
ax = nexttile(tl); hold(ax, 'on'); e = 0:7.5:90; ccol = [0.1 0.35 0.8; 0.85 0.2 0.15];
for q = 1:2, for cc = 1:2
    m = ok & C.era == er(q) & (C.MI < 0.7) == (cc == 2);
    histogram(ax, abs(C.axis_rel(m)), e, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.2, ...
        'EdgeColor', ccol(cc, :), 'LineStyle', char(ifelse(q == 1, '-', '--')));
end, end
xlabel(ax, '|preferred axis re anchor| (deg)'); ylabel(ax, 'fraction of class'); xlim(ax, [0 90]); xticks(ax, 0:15:90);
legend(ax, {'simple 2012', 'complex 2012', 'simple later', 'complex later'}, 'Box', 'off', 'Location', 'north');
title(ax, 'folded axis distribution, by rig era');
exportgraphics(f, [oa 'claude_figures/fig_types_axis_osi_mi.pdf'], 'ContentType', 'image', 'Resolution', 200); close(f);
A.n = height(C); A.cx_orth_frac = mean(abs(C.axis_rel(c)) >= 67.5); A.cx_near_frac = mean(abs(C.axis_rel(c)) <= 22.5);
A.s_near_frac = mean(abs(C.axis_rel(s)) <= 22.5); A.s_orth_frac = mean(abs(C.axis_rel(s)) >= 67.5); A.n_cx = sum(c); A.n_s = sum(s);
end
function r = ifelse(c, a, b), if c, r = a; else, r = b; end, end
