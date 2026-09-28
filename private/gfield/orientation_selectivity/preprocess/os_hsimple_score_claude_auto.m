function [sc_flank, sc_full, best_phi, Aal] = os_hsimple_score_claude_auto(T, rfimg, tmpl, loo)
% os_hsimple_score_claude_auto  Rotation-invariant match of each cell's signed RF to
% the h-simple template (tmpl from hsimple_template_2012-10-31-0.mat).
% Each RF is centered on its Vision fit (row = H - y + 0.5, col = x + 0.5; checked
% against the smoothed |z| peak: median offsets < 0.5 stixel) or, without a fit, on
% the smoothed |z| peak; resampled so its fit size (sqrt(sd1*sd2)) matches the
% template reference size (dataset median size when the cell has no fit); rotated
% through phi = 0:7.5:352.5 deg. sc_full = max_phi cosine similarity with the full
% template; sc_flank = max_phi cosine similarity with the template minus its
% radial (rotationally symmetric) average, i.e. the oriented center + flank
% structure beyond a round center. Signed: an ON-center RF scores negative against
% the OFF template. best_phi = rotation maximizing sc_flank. Aal = RF aligned at
% best_phi. loo (logical): if true, a cell that is a template member is scored
% against the template recomputed without it.
% 2026-09-28 GDF + Claude
if nargin < 4, loo = false; end
phis = 0:7.5:352.5; n = height(T);
sc_flank = nan(n, 1); sc_full = sc_flank; best_phi = sc_flank; Aal = nan(25, 25, n);
sz = sqrt(T.fit_sd1 .* T.fit_sd2); dsz = median(sz, 'omitnan');
for i = 1:n
    z = rfimg{i}; if isempty(z), continue; end
    H = size(z, 1);
    if ~isnan(T.fit_x(i)), rc = H - T.fit_y(i) + 0.5; cc = T.fit_x(i) + 0.5;
    else, zs = imgaussfilt(abs(z), 1); [~, ip] = max(zs(:)); [rc, cc] = ind2sub(size(z), ip); end
    s_i = sz(i); if isnan(s_i), s_i = dsz; end
    scl = s_i / tmpl.ref_size; if isnan(scl), scl = 1; end
    M = tmpl.M;
    if loo
        j = find(tmpl.member_ids == T.cell_id(i) & strcmp(tmpl.dataset, T.dataset{i}));
        if ~isempty(j), m = numel(tmpl.member_ids); M = (M * m - tmpl.A_members(:, :, j)) / (m - 1); end
    end
    Rm = M - radial_mean(M); M = M / norm(M(:)); Rm = Rm / norm(Rm(:));
    rf = zeros(size(phis)); rl = rf;
    for k = 1:numel(phis)
        a = os_rf_align_claude_auto(z, T.polarity(i), rc, cc, phis(k), 12, scl);
        rf(k) = a(:)' * M(:); rl(k) = a(:)' * Rm(:);
    end
    [sc_flank(i), kk] = max(rl); sc_full(i) = max(rf); best_phi(i) = phis(kk);
    Aal(:, :, i) = os_rf_align_claude_auto(z, T.polarity(i), rc, cc, phis(kk), 12, scl);
end
end
function R = radial_mean(M)
[x, y] = meshgrid(1:size(M, 2), 1:size(M, 1)); c = (size(M) + 1) / 2;
r = round(hypot(x - c(2), y - c(1))); R = zeros(size(M));
for k = unique(r(:))', R(r == k) = mean(M(r == k)); end
end
