function RF = os_rf_metrics_claude_auto(dl, gr_ids, varargin)
% os_rf_metrics_claude_auto  White-noise RF metrics for grating-run cells.
%
% Maps each grating cell to the white-noise run with map_ei (spatial EI
% correlation >= corr_thr, default 0.9, as in os_cell_finder_script) and measures
% its STA spatial RF (get_rf, color channels averaged): robust z-score
% (median/MAD over pixels); dominant polarity = sign at the |z| peak (+1 ON, -1 OFF);
% significant pixels = polarity * z > z_thr (default 4); area (pixels); number of
% connected patches (>= 2 px); fraction of significant pixels in the largest patch;
% elongation = sqrt(ratio of eigenvalues) of the z-weighted pixel covariance;
% long-axis angle (deg, STA image coordinates, 0-180). NaN when unmapped.
% 2026-09-27 GDF + Claude
p = inputParser; p.addParameter('corr_thr', 0.9); p.addParameter('z_thr', 4); p.parse(varargin{:}); o = p.Results;
gr = load_data(dl.grating_datapath); gr = load_neurons(gr); gr = load_ei(gr, gr_ids(:)');
wn = load_data(dl.wn_datapath); wn = load_neurons(wn); wn = load_ei(wn, 'all');
mp = map_ei(gr, wn, 'master_cell_type', gr_ids(:)', 'corr_threshold', o.corr_thr);
n = numel(gr_ids); z_ = nan(n, 1);
RF = table(gr_ids(:), z_, z_, z_, z_, z_, z_, z_, z_, 'VariableNames', ...
    {'cell_id', 'wn_id', 'polarity', 'area', 'n_patch', 'frac_largest', 'aspect', 'long_axis', 'peak_z'});
wid = nan(n, 1);
for i = 1:n, gi = find(gr.cell_ids == gr_ids(i)); if ~isempty(mp{gi}), wid(i) = mp{gi}; end, end
RF.wn_id = wid; ok = ~isnan(wid);
if ~any(ok), return; end
wn = load_sta(wn, 'load_sta', unique(wid(ok))');
for i = find(ok)'
    rf = get_rf(wn, wid(i)); if isempty(rf), continue; end
    rf = mean(rf, 3); z = (rf - median(rf(:))) / (1.4826 * mad(rf(:), 1));
    [pk, ip] = max(abs(z(:))); sg = sign(z(ip)); msk = sg * z > o.z_thr;
    RF.polarity(i) = sg; RF.peak_z(i) = pk; RF.area(i) = sum(msk(:));
    cc = bwconncomp(msk); sz = cellfun(@numel, cc.PixelIdxList);
    RF.n_patch(i) = sum(sz >= 2); if ~isempty(sz), RF.frac_largest(i) = max(sz) / sum(sz); end
    [yy, xx] = find(msk); if numel(xx) < 3, continue; end
    w = sg * z(msk); w = w / sum(w); Xc = [xx yy] - [sum(w .* xx) sum(w .* yy)];
    C = (Xc .* w)' * Xc; [V, E] = eig(C); [ev, io] = sort(diag(E), 'descend');
    RF.aspect(i) = sqrt(ev(1) / max(ev(2), eps)); RF.long_axis(i) = mod(atan2d(V(2, io(1)), V(1, io(1))), 180);
end
RF.Properties.UserData = struct('rf_size', size(get_rf(wn, wid(find(ok, 1)))), 'wn_path', dl.wn_datapath);
end
