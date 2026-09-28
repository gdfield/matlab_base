function out = os_types_hsimple_claude_auto(dataset, varargin)
% os_types_hsimple_claude_auto  h-simple anchored OS typing for one dataset (v3).
%
% Cells: curated OS cells (old: GDF definite + maybe; new: provisional classifier
% definite + maybe) from os_cells_master_table.mat.
%  1. RF polarity from the STA time course (significant_stixels ->
%     time_course_from_sta -> guess_polarity), Vision params time course when that
%     is 0/NaN. Signed RF = polarity * robust z of get_rf.
%  2. h-simple template score (os_hsimple_score_claude_auto): rotation-invariant
%     match of the signed, size-normalized RF to the 2012-10-31-0 h-simple template
%     (flank score = oriented center + flank structure beyond a round center).
%  3. Simple = MI >= mi_thr. Candidate modes = local maxima of an axial von Mises
%     kernel density (kappa on 2*theta) of simple-cell preferred axes with >= 3 cells
%     within +/- win. Anchor = mode whose cells have the highest median flank score
%     (>= 3 cells with RFs); if no mode qualifies, the mode with most cells
%     (method 'count (no RFs)'). Margin and one-sided rank-sum p vs the next mode kept.
%  4. h-simple = simple cells within +/- win of the anchor. Complex vOS = cells within
%     +/- win of anchor + 90 with MI < mi_thr (v3: MI only; GDF 2026-09-28, RF patch
%     rule dropped as unreliable). Orth-simple = orthogonal cells with MI >= mi_thr.
%  5. Dataset excluded from further consideration if n(h-simple) < min_h (5).
% v1/v2 (2026-09-27) chose the anchor by cell count and used an RF-patch rule.
% Writes os_analysis/claude_types/<dataset>/ : cells.csv, summary.pdf,
% mosaics.pdf, rf_gallery_<class>.pdf, tuning_polar.pdf.
% 2026-09-27/28 GDF + Claude
p = inputParser;
p.addParameter('mi_thr', 0.7); p.addParameter('win', 22.5); p.addParameter('kappa', 8);
p.addParameter('amb_ratio', 0.6); p.addParameter('patch_frac', 0.6); p.addParameter('sd_scale', 1);
p.addParameter('do_figs', true); p.addParameter('min_h', 5);
p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
od = [oa 'claude_types/' dataset '/']; if ~exist(od, 'dir'), mkdir(od); end
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS');
T = S.OS(strcmp(S.OS.dataset, dataset), :);
reg = load([oa 'claude_features/registry_used.mat']); dl = reg.allds(strcmp(reg.nm, dataset));
rff = [oa 'claude_rf/' dataset '_rf.mat']; has_rf = exist(rff, 'file') == 2;
T.wn_id = nan(height(T), 1);
if has_rf, has_rf = ismember('RF', who('-file', rff)); end
if has_rf, R = load(rff, 'RF'); [~, loc] = ismember(T.cell_id, R.RF.cell_id); T.wn_id(loc > 0) = R.RF.wn_id(loc(loc > 0)); end
has_rf = has_rf && any(~isnan(T.rf_area));
ad = @(a, b) abs(mod(a - b + 90, 180) - 90);   % axial difference, 0-90
% ---- 4 polarity, RF images, fits
T.polarity = nan(height(T), 1); T.polarity_vision = nan(height(T), 1);
T.fit_x = nan(height(T), 1); T.fit_y = T.fit_x; T.fit_sd1 = T.fit_x; T.fit_sd2 = T.fit_x; T.fit_ang = T.fit_x;
rfimg = cell(height(T), 1); wn = [];
ok = ~isnan(T.wn_id);
if any(ok)
    wn = load_data(dl.wn_datapath); wn = load_neurons(wn);
    try, wn = load_params(wn); catch, end
    wn = load_sta(wn, 'load_sta', unique(T.wn_id(ok))');
    for i = find(ok)'
        idx = get_cell_indices(wn, T.wn_id(i)); sta = wn.stas.stas{idx};
        try
            sg = significant_stixels(sta); tcv = time_course_from_sta(sta, sg); T.polarity(i) = guess_polarity(tcv);
        catch, end
        if isfield(wn, 'vision') && isfield(wn.vision, 'timecourses')
            v = wn.vision.timecourses(idx); tv = [v.r(:) v.g(:) v.b(:)]; T.polarity_vision(i) = guess_polarity(tv);
        end
        if isfield(wn, 'vision') && isfield(wn.vision, 'sta_fits') && ~isempty(wn.vision.sta_fits{idx})
            f = wn.vision.sta_fits{idx}; T.fit_x(i) = f.mean(1); T.fit_y(i) = f.mean(2);
            T.fit_sd1(i) = f.sd(1); T.fit_sd2(i) = f.sd(2); T.fit_ang(i) = f.angle;
        end
        rf = get_rf(wn, T.wn_id(i));
        if ~isempty(rf), rf = mean(rf, 3); rfimg{i} = (rf - median(rf(:))) / (1.4826 * mad(rf(:), 1)); end
    end
end
% polarity: STA time course; if 0/NaN, use Vision params time course (v2)
T.polarity_source = repmat("sta", height(T), 1);
fb = (isnan(T.polarity) | T.polarity == 0) & ~isnan(T.polarity_vision) & T.polarity_vision ~= 0;
T.polarity(fb) = T.polarity_vision(fb); T.polarity_source(fb) = "vision"; T.polarity_source(isnan(T.polarity) | T.polarity == 0) = "none";
% ---- h-simple template score (v3)
tp = load([oa 'claude_types/hsimple_template_2012-10-31-0.mat'], 'tmpl');
[T.flank_score, T.full_score, T.best_phi] = os_hsimple_score_claude_auto(T, rfimg, tp.tmpl, true);
% ---- anchor (v3): candidate modes of simple-cell preferred axes; pick the mode whose
% cells best match the h-simple template (median flank score, >= 3 cells with RFs)
ad = @(a, b) abs(mod(a - b + 90, 180) - 90);   % axial difference, 0-90
simple = T.MI >= o.mi_thr; elig = simple & ~isnan(T.pref_axis);
th = 0:0.5:179.5; kd = zeros(size(th));
for k = find(elig)', kd = kd + exp(o.kappa * cosd(2 * (th - T.pref_axis(k)))); end
pk = th(kd >= circshift(kd, 1) & kd >= circshift(kd, -1) & kd > 0);
M_ = struct('axis', {}, 'n', {}, 'n_rf', {}, 'med_score', {}, 'scores', {});
for q = pk
    m = elig & ad(T.pref_axis, q) <= o.win; if sum(m) < 3, continue; end
    sc = T.flank_score(m & ~isnan(T.flank_score));
    M_(end + 1) = struct('axis', q, 'n', sum(m), 'n_rf', numel(sc), 'med_score', median(sc), 'scores', sc); %#ok<AGROW>
end
anc = NaN; method = "none"; margin = NaN; p_sep = NaN; second = NaN;
if ~isempty(M_)
    okm = [M_.n_rf] >= 3;
    if any(okm)
        ms = [M_.med_score]; ms(~okm) = -inf; [~, b] = max(ms); anc = M_(b).axis; method = "template";
        ms(b) = -inf; if any(isfinite(ms)), [~, b2] = max(ms); second = M_(b2).axis;
            margin = M_(b).med_score - M_(b2).med_score; p_sep = ranksum(M_(b).scores, M_(b2).scores, 'tail', 'right'); end
    else
        [~, b] = max([M_.n]); anc = M_(b).axis; method = "count (no RFs)";
    end
end
orth = mod(anc + 90, 180);
nearA = ad(T.pref_axis, anc) <= o.win; nearO = ad(T.pref_axis, orth) <= o.win;
cls = repmat("other", height(T), 1); reason = repmat("", height(T), 1);
cls(elig & nearA) = "h_simple";
lowMI = T.MI < o.mi_thr;
cv = nearO & lowMI; cls(cv) = "complex_vOS"; reason(cv) = "MI";
cls(nearO & ~cv & simple) = "orth_simple";
nA = sum(cls == "h_simple"); nOe = sum(elig & nearO);
amb = nOe >= o.amb_ratio * max(nA, 1);
T.class = cls; T.reason = reason; T.axis_rel_anchor = mod(T.pref_axis - anc + 90, 180) - 90;
% EI fallback positions (grating run) for cells without an RF fit
T.ei_x = nan(height(T), 1); T.ei_y = T.ei_x;
if any(isnan(T.fit_x))
    try
        gr = load_data(dl.grating_datapath); gr = load_neurons(gr); gr = load_ei(gr, T.cell_id(:)');
        pos = gr.ei.position;
        for i = 1:height(T)
            try, ei = gr.ei.eis{get_cell_indices(gr, T.cell_id(i))}; catch, continue; end
            if isempty(ei), continue; end
            a = max(ei, [], 2) - min(ei, [], 2); [~, s] = sort(a, 'descend'); s = s(1:min(5, numel(s)));
            w = a(s) .^ 2; T.ei_x(i) = sum(w .* pos(s, 1)) / sum(w); T.ei_y(i) = sum(w .* pos(s, 2)) / sum(w);
        end
    catch ME, warning('EI fallback failed for %s: %s', dataset, ME.message); end
end
% ---- tuning curves (best condition)
Fz = load([oa 'claude_features/' dataset '_features.mat'], 'F'); Fz = Fz.F;
T.tuning = nan(height(T), numel(Fz.dirs));
for i = 1:height(T)
    r = find(Fz.cell_ids == T.cell_id(i), 1); is = find(Fz.sps == T.best_sp(i), 1); it = find(Fz.tps == T.best_tp(i), 1);
    if ~isempty(r) && ~isempty(is) && ~isempty(it), T.tuning(i, :) = Fz.tuning{is, it}(r, :); end
end
Tout = removevars(T, 'tuning'); writetable(Tout, [od 'cells.csv']);
info = struct('dataset', dataset, 'anchor_axis', anc, 'anchor_method', method, 'anchor_second_mode', second, ...
    'anchor_margin', margin, 'anchor_p_sep', p_sep, 'modes', rmfield(M_, 'scores'), 'excluded', nA < o.min_h, 'n_h_simple', nA, 'n_complex_vOS', sum(cls == "complex_vOS"), ...
    'n_orth_simple', sum(cls == "orth_simple"), 'n_other', sum(cls == "other"), 'n_orth_eligible', nOe, ...
    'ambiguous_anchor', amb, 'has_rf', has_rf, 'n_cells', height(T), 'params', o, 'version', 3);
save([od 'types.mat'], 'T', 'info', 'rfimg');
out = info;
if o.do_figs, os_types_figs_claude_auto(T, info, rfimg, Fz.dirs, od, o); end
end
