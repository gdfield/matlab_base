function os_types_figs_claude_auto(T, info, rfimg, dirs, od, o)
% os_types_figs_claude_auto  Figures for os_types_hsimple_claude_auto (one dataset).
% summary.pdf: preferred-axis distributions and MI vs axis relative to the anchor.
% mosaics.pdf: RF fits (sd_scale x SD ellipses, Vision coordinates, stixels, y up)
%   of h-simple and complex vOS cells; other OS cells grey; ON solid, OFF dashed.
%   Cells without a fit shown as EI centers (grating-run electrode coordinates, um)
%   in a separate panel.
% rf_gallery_<class>.pdf: signed RFs (polarity x robust z; white = ON, black = OFF),
%   25x25-stixel crop around the |z| peak, preferred drift axis in orange drawn
%   with y up (v3; grating and STA coordinates assumed to match), sorted by
%   h-simple template flank score. summary.pdf panel 3 (v3): flank score vs axis.
% tuning_polar.pdf: normalized best-condition direction tuning, by class.
% 2026-09-27 GDF + Claude
ds = info.dataset; C = ["h_simple", "complex_vOS", "orth_simple"];
col = struct('h_simple', [0.1 0.35 0.8], 'complex_vOS', [0.85 0.2 0.15], 'orth_simple', [0.2 0.6 0.3], 'other', [0.6 0.6 0.6]);
ttl = sprintf('%s | anchor %.0f deg (%s, margin %.2f, p %.2g) | h-simple %d, complex vOS %d, orth-simple %d, other %d%s', ds, info.anchor_axis, ...
    info.anchor_method, info.anchor_margin, info.anchor_p_sep, info.n_h_simple, info.n_complex_vOS, info.n_orth_simple, info.n_other, ...
    string(ifelse(info.excluded, ' | EXCLUDED (<5 h-simple)', '')));
% ---------- summary
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 3.6], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); e = 0:15:180; hold(ax, 'on');
s = T.MI >= o.mi_thr; histogram(ax, T.pref_axis(s), e, 'FaceColor', col.h_simple, 'FaceAlpha', 0.5);
histogram(ax, T.pref_axis(~s & ~isnan(T.MI)), e, 'FaceColor', col.complex_vOS, 'FaceAlpha', 0.5);
yl = ylim(ax); xline(ax, info.anchor_axis, 'k-', 'LineWidth', 1.2); xline(ax, mod(info.anchor_axis + 90, 180), 'k--', 'LineWidth', 1.2);
xlim(ax, [0 180]); xlabel(ax, 'preferred drift axis (deg, stimulus coords)'); ylabel(ax, 'cells');
legend(ax, {sprintf('MI>=%.1f', o.mi_thr), sprintf('MI<%.1f', o.mi_thr), 'anchor', 'anchor+90'}, 'Box', 'off', 'Location', 'best'); title(ax, 'all curated OS cells');
ax = nexttile(tl); hold(ax, 'on');
for c = ["other", C], m = T.class == c; scatter(ax, T.axis_rel_anchor(m), T.MI(m), 14, col.(c), 'filled', 'MarkerFaceAlpha', 0.7); end
xline(ax, [-o.win o.win], 'k:'); xline(ax, [90 - o.win, -90 + o.win], 'k:'); yline(ax, o.mi_thr, 'k:');
xlim(ax, [-90 90]); xlabel(ax, 'preferred axis - anchor (deg)'); ylabel(ax, 'MI = 2F1/F0'); title(ax, 'class assignment');
legend(ax, {'other', 'h-simple', 'complex vOS', 'orth-simple'}, 'Box', 'off', 'Location', 'best');
ax = nexttile(tl); hold(ax, 'on');
for c = ["other", C], m = T.class == c; scatter(ax, T.axis_rel_anchor(m), T.flank_score(m), 14, col.(c), 'filled', 'MarkerFaceAlpha', 0.7); end
xline(ax, [-o.win o.win], 'k:'); xlim(ax, [-90 90]); xlabel(ax, 'preferred axis - anchor (deg)'); ylabel(ax, 'h-simple template flank score');
title(ax, 'RF match to h-simple template');
title(tl, ttl, 'FontSize', 8, 'Interpreter', 'none');
exportgraphics(f, [od 'summary.pdf'], 'ContentType', 'vector'); close(f);
% ---------- mosaics
hasfit = ~isnan(T.fit_x); hasei = ~hasfit & ~isnan(T.ei_x);
np = 2 + any(hasei);
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 4 * np 4], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, np, 'TileSpacing', 'compact', 'Padding', 'compact'); tt = linspace(0, 2 * pi, 60);
for c = C(1:2)
    ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
    for i = find(hasfit)'
        a = T.fit_ang(i); sx = o.sd_scale * T.fit_sd1(i); sy = o.sd_scale * T.fit_sd2(i);
        x = T.fit_x(i) + sx * cos(tt) * cos(a) - sy * sin(tt) * sin(a); y = T.fit_y(i) + sx * cos(tt) * sin(a) + sy * sin(tt) * cos(a);
        if T.class(i) == c
            ls = '-'; if T.polarity(i) < 0, ls = '--'; elseif ~(T.polarity(i) > 0), ls = ':'; end
            plot(ax, x, y, ls, 'Color', col.(c), 'LineWidth', 1.1);
        else, plot(ax, x, y, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4); end
    end
    m = hasfit & T.class == c;
    title(ax, sprintf('%s (n=%d with fit; ON %d, OFF %d)', strrep(c, '_', ' '), sum(m), sum(m & T.polarity > 0), sum(m & T.polarity < 0)));
    xlabel(ax, 'stixels'); ylabel(ax, 'stixels'); box(ax, 'on');
end
if any(hasei)
    ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
    for c = ["other", C], m = hasei & T.class == c; scatter(ax, T.ei_x(m), T.ei_y(m), 16, col.(c), 'filled'); end
    title(ax, sprintf('cells without RF fit: EI centers (n=%d); blue h-simple, red complex vOS, green orth-simple, grey other', sum(hasei))); xlabel(ax, 'um'); ylabel(ax, 'um'); box(ax, 'on');
end
title(tl, [ttl sprintf(' | %g-SD ellipses; solid ON, dashed OFF, dotted unknown; grey = other OS cells', o.sd_scale)], 'FontSize', 7, 'Interpreter', 'none');
exportgraphics(f, [od 'mosaics.pdf'], 'ContentType', 'vector'); close(f);
% ---------- RF galleries
for c = C
    idx = find(T.class == c & ~cellfun(@isempty, rfimg)); fn = [od 'rf_gallery_' char(c) '.pdf'];
    if exist(fn, 'file'), delete(fn); end
    if isempty(idx), continue; end
    n_norf = sum(T.class == c) - numel(idx);
    [~, so] = sort(T.flank_score(idx), 'descend'); idx = idx(so); per = 48; np_ = ceil(numel(idx) / per);
    for pg = 1:np_
        f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 11 8.5], 'DefaultAxesFontSize', 6);
        tl = tiledlayout(f, 6, 8, 'TileSpacing', 'tight', 'Padding', 'compact');
        for k = (pg - 1) * per + 1:min(pg * per, numel(idx))
            i = idx(k); z = rfimg{i}; pol = T.polarity(i); if isnan(pol) || pol == 0, pol = 1; end
            zs = pol * z; [~, ip] = max(abs(z(:))); [yy, xx] = ind2sub(size(z), ip); hw = 12;
            r1 = max(1, yy - hw):min(size(z, 1), yy + hw); c1 = max(1, xx - hw):min(size(z, 2), xx + hw); zc = zs(r1, c1);
            ax = nexttile(tl); imagesc(ax, zc); axis(ax, 'image', 'off'); colormap(ax, gray); mx = max(abs(zc(:)));
            clim(ax, [-mx mx]); hold(ax, 'on'); yy = yy - r1(1) + 1; xx = xx - c1(1) + 1; L_ = 6;
            plot(ax, xx + L_ * cosd(T.pref_axis(i)) * [-1 1], yy - L_ * sind(T.pref_axis(i)) * [-1 1], '-', 'Color', [0.95 0.5 0.1], 'LineWidth', 1.2);
            pstr = 'ON'; if T.polarity(i) < 0, pstr = 'OFF'; elseif isnan(T.polarity(i)) || T.polarity(i) == 0, pstr = '?'; end
            xlim(ax, [0.5 size(zc, 2) + 0.5]); ylim(ax, [0.5 size(zc, 1) + 0.5]);
            title(ax, sprintf('%d %s MI%.2f\\newlineax%.0f fs%.2f', T.cell_id(i), pstr, T.MI(i), T.pref_axis(i), T.flank_score(i)), 'FontSize', 5, 'FontWeight', 'normal');
        end
        title(tl, sprintf('%s | %s | page %d/%d | %d cells shown, %d without RF not shown | signed RF (white ON, black OFF), up to 25x25 stixels, orange = preferred drift axis (y up); sorted by template flank score (fs)', ds, strrep(c, '_', ' '), pg, np_, numel(idx), n_norf), 'FontSize', 8, 'Interpreter', 'none');
        exportgraphics(f, fn, 'ContentType', 'image', 'Resolution', 150, 'Append', pg > 1); close(f);
    end
end
% ---------- polar tuning
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 9 3.2], 'DefaultAxesFontSize', 7);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); th = deg2rad([dirs(:); dirs(1)]);
for c = C
    pax = polaraxes(tl); pax.Layout.Tile = find(C == c); hold(pax, 'on'); m = find(T.class == c);
    for i = m', r = T.tuning(i, :); if all(isnan(r)) || max(r) <= 0, continue; end
        r = r / max(r); polarplot(pax, th, [r(:); r(1)], '-', 'Color', [col.(c) 0.35], 'LineWidth', 0.6); end
    a = deg2rad(info.anchor_axis); polarplot(pax, [a a + pi], [1 1], 'k-', 'LineWidth', 1.2);
    polarplot(pax, [a + pi / 2, a + 3 * pi / 2], [1 1], 'k--', 'LineWidth', 1.2); rlim(pax, [0 1]);
    title(pax, sprintf('%s (n=%d)', strrep(c, '_', ' '), numel(m)));
end
title(tl, [ttl ' | best condition, normalized; solid = anchor axis, dashed = orthogonal'], 'FontSize', 7, 'Interpreter', 'none');
exportgraphics(f, [od 'tuning_polar.pdf'], 'ContentType', 'vector'); close(f);
end
function r = ifelse(c, a, b), if c, r = a; else, r = b; end, end
