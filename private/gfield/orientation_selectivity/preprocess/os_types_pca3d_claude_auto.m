function [fig, P] = os_types_pca3d_claude_auto(varargin)
% os_types_pca3d_claude_auto  PCA of four per-cell features for curated OS cells in
% the 36 kept M9 datasets, shown interactively in the space of PCs 1-3.
% Features (each z-scored across cells):
%   axis_c  = cos(2 * axis re h-simple anchor): +1 at the anchor, -1 orthogonal
%             (orientation is axial, so the angle enters through 2*theta; the sign of
%             the angle, which a mirror could flip between rigs, drops out)
%   OSI     = pipeline osi_max (max 0.5)
%   MI      = 2F1/F0
%   rel     = rel_pref from os_pref_reliability_claude_auto (mean pairwise correlation
%             of 50-ms PSTHs across repetitions, 1-8 s, preferred direction, best
%             condition). Option 'rel_measure', 'fano' uses -log10(Fano factor of the
%             preferred-direction count) instead (count reliability, not temporal).
% Colors: blue/light blue MI >= 0.7 near anchor / near orthogonal; red/orange MI < 0.7
% near anchor / near orthogonal. Circles 2012 rig, triangles later rig. Black lines =
% feature loadings (biplot, scaled). Data tips give dataset and cell id.
% Saves claude_figures/fig_types_pca3d[_fano].fig/.pdf and claude_types/pca_features[_fano].csv.
% 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('rel_measure', 'psth'); p.addParameter('visible', 'on'); p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
C = readtable([oa 'claude_types/axis_osi_mi_cells.csv'], 'TextType', 'string'); C.dataset = string(C.dataset); C.class = string(C.class);
dd = dir([oa 'claude_types/pref_reliability/*.mat']); Rs = cell(numel(dd), 1);
for k = 1:numel(dd), q = load(fullfile(dd(k).folder, dd(k).name), 'Rl'); Rs{k} = q.Rl; end
C = innerjoin(C, vertcat(Rs{:}), 'Keys', {'dataset', 'cell_id'});
C.era = repmat("later", height(C), 1); C.era(startsWith(C.dataset, "2012")) = "2012";
if strcmp(o.rel_measure, 'fano'), rel = -log10(C.fano_pref); rn = '-log10 Fano'; sfx = '_fano';
else, rel = C.rel_pref; rn = 'PSTH reliability'; sfx = ''; end
X = [cosd(2 * C.axis_rel), C.OSI, C.MI, rel]; fn = {'cos 2*axis', 'OSI', 'MI', rn};
ok = all(isfinite(X), 2); X = X(ok, :); C = C(ok, :);
Z = (X - mean(X)) ./ std(X);
[coeff, score, ~, ~, expl] = pca(Z);
C.PC1 = score(:, 1); C.PC2 = score(:, 2); C.PC3 = score(:, 3); C.rel_used = X(:, 4);
writetable(C, [oa 'claude_types/pca_features' sfx '.csv']);
nearO = abs(C.axis_rel) > 45;
G = {C.MI >= 0.7 & ~nearO, C.MI >= 0.7 & nearO, C.MI < 0.7 & ~nearO, C.MI < 0.7 & nearO};
gc = [0.1 0.35 0.8; 0.45 0.75 1; 0.85 0.2 0.15; 1 0.6 0.2];
gn = {'MI>=0.7 near anchor', 'MI>=0.7 near orthogonal', 'MI<0.7 near anchor', 'MI<0.7 near orthogonal'};
mk = {'o', '^'}; er = ["2012", "later"];
fig = figure('Name', ['PCA of axis, OSI, MI, ' rn], 'Color', 'w', 'Visible', o.visible);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); nm = {};
for e = 1:2
    for c = 1:4
        m = C.era == er(e) & G{c}; if ~any(m), continue; end
        h = scatter3(ax, C.PC1(m), C.PC2(m), C.PC3(m), 20, gc(c, :), mk{e}, 'filled', 'MarkerFaceAlpha', 0.6);
        h.DataTipTemplate.DataTipRows = [dataTipTextRow('dataset', C.dataset(m)), dataTipTextRow('cell', C.cell_id(m)), ...
            dataTipTextRow('class', C.class(m)), dataTipTextRow('axis re anchor', C.axis_rel(m)), dataTipTextRow('MI', C.MI(m)), ...
            dataTipTextRow('OSI', C.OSI(m)), dataTipTextRow(rn, C.rel_used(m))];
        nm{end + 1} = sprintf('%s, %s (n=%d)', gn{c}, er(e), sum(m)); %#ok<AGROW>
    end
end
sc = 0.9 * max(abs(score(:, 1:3)), [], 'all');
for f = 1:4
    plot3(ax, [0 sc * coeff(f, 1)], [0 sc * coeff(f, 2)], [0 sc * coeff(f, 3)], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(ax, 1.08 * sc * coeff(f, 1), 1.08 * sc * coeff(f, 2), 1.08 * sc * coeff(f, 3), fn{f}, 'FontWeight', 'bold');
end
xlabel(ax, sprintf('PC1 (%.0f%%)', expl(1))); ylabel(ax, sprintf('PC2 (%.0f%%)', expl(2))); zlabel(ax, sprintf('PC3 (%.0f%%)', expl(3)));
legend(ax, nm, 'Location', 'northeastoutside', 'Box', 'off');
title(ax, sprintf('%d cells, 36 datasets; features z-scored; drag to rotate', height(C)));
view(ax, -35, 25); rotate3d(fig, 'on');
savefig(fig, [oa 'claude_figures/fig_types_pca3d' sfx '.fig']);
exportgraphics(fig, [oa 'claude_figures/fig_types_pca3d' sfx '.pdf'], 'ContentType', 'image', 'Resolution', 200);
P = struct('coeff', coeff, 'explained', expl, 'features', {fn}, 'n', height(C), 'corr', corr(X));
end
