function os_rf_gallery_claude_auto(dataset, n_each)
% os_rf_gallery_claude_auto  White-noise RFs of the most orientation-tuned simple
% (MI >= 0.7) and complex (MI < 0.7) OS cells in one dataset, with each cell's
% preferred drift axis (stimulus coordinates) drawn for reference. Writes
% os_analysis/claude_figures/fig_rf_gallery_<dataset>.pdf.
% 2026-09-27 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); OS = S.OS;
R = load([oa 'claude_rf/' dataset '_rf.mat'], 'RF'); RF = R.RF;
reg = load([oa 'claude_features/registry_used.mat']); dl = reg.allds(strcmp(reg.nm, dataset));
s = strcmp(OS.dataset, dataset) & ~isnan(OS.rf_area);
T = OS(s, :); [~, loc] = ismember(T.cell_id, RF.cell_id); T.wn_id = RF.wn_id(loc);
Ts = sortrows(T(~T.cx, :), 'F2n', 'descend'); Tc = sortrows(T(T.cx, :), 'F2n', 'descend');
Ts = Ts(1:min(n_each, height(Ts)), :); Tc = Tc(1:min(n_each, height(Tc)), :);
wn = load_data(dl.wn_datapath); wn = load_neurons(wn); wn = load_sta(wn, 'load_sta', [Ts.wn_id; Tc.wn_id]');
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 3.2], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 2, n_each, 'TileSpacing', 'tight', 'Padding', 'compact');
rows = {Ts, Tc}; nmr = {'simple', 'complex'};
for r = 1:2
    for i = 1:n_each
        ax = nexttile(tl, (r - 1) * n_each + i); axis(ax, 'off');
        if i > height(rows{r}), continue; end
        c = rows{r}(i, :); rf = mean(get_rf(wn, c.wn_id), 3); z = (rf - median(rf(:))) / (1.4826 * mad(rf(:), 1));
        [~, ip] = max(abs(z(:))); [yy, xx] = ind2sub(size(z), ip); hw = 12;
        r1 = max(1, yy - hw):min(size(z, 1), yy + hw); c1 = max(1, xx - hw):min(size(z, 2), xx + hw); zc = z(r1, c1);
        imagesc(ax, zc); axis(ax, 'image', 'off'); colormap(ax, gray); hold(ax, 'on'); clim(ax, [-max(abs(zc(:))) max(abs(zc(:)))]);
        yy = yy - r1(1) + 1; xx = xx - c1(1) + 1; L_ = 6;
        plot(ax, xx + L_ * cosd(c.pref_axis) * [-1 1], yy + L_ * sind(c.pref_axis) * [-1 1], '-', 'Color', [0.9 0.3 0.1], 'LineWidth', 1.5);
        title(ax, sprintf('%s %d MI %.2f\\newlineF2n %.2f axis %.0f', nmr{r}, c.cell_id, c.MI, c.F2n, c.pref_axis), 'FontSize', 6, 'FontWeight', 'normal');
    end
end
title(tl, sprintf('%s: top row simple (MI>=0.7), bottom row complex (MI<0.7); 25x25-stixel crop around RF peak; orange = preferred drift axis', dataset), 'FontSize', 8);
exportgraphics(f, [oa 'claude_figures/fig_rf_gallery_' dataset '.pdf'], 'ContentType', 'image', 'Resolution', 200); close(f);
end
