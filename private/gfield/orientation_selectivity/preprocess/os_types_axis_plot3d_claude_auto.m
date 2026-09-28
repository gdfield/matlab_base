function fig = os_types_axis_plot3d_claude_auto(varargin)
% os_types_axis_plot3d_claude_auto  Interactive 3D scatter of preferred axis
% (relative to each dataset's h-simple anchor, M9 v3), MI and OSI for curated OS
% cells in the 36 kept datasets (data: claude_types/axis_osi_mi_cells.csv, written
% by os_types_axis_plot_claude_auto). Opens a visible figure with rotate3d on;
% data tips show dataset, cell id, class and values. Optional 'fold' (false):
% plot |relative axis| (0-90) instead of -90..90. Optional 'wrap45' (false): wrap the
% relative axis modulo 90 into -45..45, so anchor-aligned and orthogonal cells both sit
% at 0; colors then also mark which side (near anchor / near orthogonal) each came from. Color: blue MI >= 0.7 (simple),
% red MI < 0.7 (complex); circles = 2012 rig, triangles = later rig.
% Also saves claude_figures/fig_types_axis_osi_mi_3d.fig (reopen with openfig).
% 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('fold', false); p.addParameter('wrap45', false); p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
C = readtable([oa 'claude_types/axis_osi_mi_cells.csv'], 'TextType', 'string');
C = C(~isnan(C.axis_rel) & ~isnan(C.MI) & ~isnan(C.OSI), :);
C.dataset = string(C.dataset); C.class = string(C.class);
C.era = repmat("later", height(C), 1); C.era(startsWith(C.dataset, "2012")) = "2012";   % era column re-derived (readtable parses it as numeric)
x = C.axis_rel; if o.fold, x = abs(x); end
nearO = abs(C.axis_rel) > 45;   % closer to the orthogonal axis than to the anchor
if o.wrap45, x = mod(C.axis_rel + 45, 90) - 45; end   % anchor and orthogonal axes both map to 0
fig = figure('Name', 'OSI vs MI vs axis (re h-simple anchor)', 'Color', 'w', 'Visible', 'on');
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
cl = [0.1 0.35 0.8; 0.85 0.2 0.15]; mk = {'o', '^'}; er = ["2012", "later"]; nm = {};
if o.wrap45   % groups: simple/complex x near anchor/near orthogonal; markers by rig
    G = {C.MI >= 0.7 & ~nearO, C.MI >= 0.7 & nearO, C.MI < 0.7 & ~nearO, C.MI < 0.7 & nearO};
    gc = [0.1 0.35 0.8; 0.45 0.75 1; 0.85 0.2 0.15; 1 0.6 0.2];
    gn = {'MI>=0.7 near anchor', 'MI>=0.7 near orthogonal', 'MI<0.7 near anchor', 'MI<0.7 near orthogonal'};
else
    G = {C.MI >= 0.7, C.MI < 0.7}; gc = cl; gn = {'MI>=0.7', 'MI<0.7'};
end
for e = 1:2
    for c = 1:numel(G)
        m = C.era == er(e) & G{c}; if ~any(m), continue; end
        h = scatter3(ax, x(m), C.MI(m), C.OSI(m), 22, gc(c, :), mk{e}, 'filled', 'MarkerFaceAlpha', 0.6);
        h.DataTipTemplate.DataTipRows = [dataTipTextRow('dataset', C.dataset(m)), dataTipTextRow('cell', C.cell_id(m)), ...
            dataTipTextRow('class', C.class(m)), dataTipTextRow('axis re anchor', C.axis_rel(m)), ...
            dataTipTextRow('MI', C.MI(m)), dataTipTextRow('OSI', C.OSI(m))];
        nm{end + 1} = sprintf('%s, %s (n=%d)', gn{c}, er(e), sum(m)); %#ok<AGROW>
    end
end
if o.wrap45, xlim(ax, [-45 45]); xticks(ax, -45:15:45); xlabel(ax, 'preferred axis re anchor, wrapped mod 90 (deg; anchor and orthogonal both at 0)');
elseif o.fold, xlim(ax, [0 90]); xticks(ax, 0:15:90); xlabel(ax, '|preferred axis re h-simple anchor| (deg)');
else, xlim(ax, [-90 90]); xticks(ax, -90:45:90); xlabel(ax, 'preferred axis re h-simple anchor (deg)'); end
ylabel(ax, 'MI (2F1/F0)'); zlabel(ax, 'OSI (pipeline, max 0.5)');
legend(ax, nm, 'Location', 'northeastoutside', 'Box', 'off');
title(ax, sprintf('%d curated OS cells, %d datasets; drag to rotate', height(C), numel(unique(C.dataset))));
view(ax, -35, 25); rotate3d(fig, 'on');
sfx = ''; if o.fold, sfx = '_folded'; end; if o.wrap45, sfx = '_wrap45'; end
savefig(fig, [oa 'claude_figures/fig_types_axis_osi_mi_3d' sfx '.fig']);
end
function r = ifelse(c, a, b), if c, r = a; else, r = b; end, end
