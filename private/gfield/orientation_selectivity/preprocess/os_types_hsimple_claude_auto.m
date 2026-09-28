function out = os_types_hsimple_claude_auto(dataset, varargin)
% os_types_hsimple_claude_auto  h-simple anchored OS typing for one dataset.
%
% Cells: curated OS cells (old: GDF definite + maybe; new: provisional classifier
% definite + maybe) from os_cells_master_table.mat.
%  1. Anchor-eligible = MI >= mi_thr (simple) AND, when a WN RF exists, a single
%     RF patch (rf_npatch == 1). Cells without an RF (or datasets without WN RFs)
%     are judged on MI alone (v2, after QC: v1 excluded them when the dataset had RFs).
%  2. Anchor axis = peak of an axial von Mises kernel density (kappa on 2*theta)
%     of anchor-eligible preferred axes. h-simple = eligible cells within
%     +/- win deg of the anchor axis. Ambiguity flag if the eligible cells
%     within +/- win of the orthogonal axis number >= amb_ratio * n(h-simple).
%  3. Complex vOS candidate = preferred axis within +/- win of anchor + 90 AND
%     (MI < mi_thr OR rf_npatch >= 2 OR rf_frac_largest < patch_frac).
%     Orthogonal cells failing both = 'orth_simple'. All else = 'other'.
%  4. RF polarity from the STA time course (significant_stixels ->
%     time_course_from_sta -> guess_polarity); Vision params time course polarity
%     kept as a cross-check. Signed RF = polarity * robust z of get_rf.
% Writes os_analysis/claude_types/<dataset>/ : cells.csv, summary.pdf,
% mosaics.pdf, rf_gallery_<class>.pdf, tuning_polar.pdf.
% 2026-09-27 GDF + Claude
p = inputParser;
p.addParameter('mi_thr', 0.7); p.addParameter('win', 22.5); p.addParameter('kappa', 8);
p.addParameter('amb_ratio', 0.6); p.addParameter('patch_frac', 0.6); p.addParameter('sd_scale', 1);
p.addParameter('do_figs', true);
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
% ---- 1-2 anchor
simple = T.MI >= o.mi_thr;
elig = simple & ~isnan(T.pref_axis) & (isnan(T.rf_npatch) | T.rf_npatch == 1);
th = 0:0.5:179.5; kd = zeros(size(th));
for k = find(elig)', kd = kd + exp(o.kappa * cosd(2 * (th - T.pref_axis(k)))); end
[~, im] = max(kd); anc = th(im); if ~any(elig), anc = NaN; end
orth = mod(anc + 90, 180);
nearA = ad(T.pref_axis, anc) <= o.win; nearO = ad(T.pref_axis, orth) <= o.win;
cls = repmat("other", height(T), 1); reason = repmat("", height(T), 1);
cls(elig & nearA) = "h_simple";
lowMI = T.MI < o.mi_thr; patchy = T.rf_npatch >= 2 | T.rf_frac_largest < o.patch_frac;
cv = nearO & (lowMI | patchy); cls(cv) = "complex_vOS";
reason(cv & lowMI & patchy) = "MI+RF"; reason(cv & lowMI & ~patchy) = "MI"; reason(cv & ~lowMI & patchy) = "RF";
cls(nearO & ~cv) = "orth_simple";
nA = sum(cls == "h_simple"); nOe = sum(elig & nearO);
amb = nOe >= o.amb_ratio * max(nA, 1);
T.class = cls; T.reason = reason; T.axis_rel_anchor = mod(T.pref_axis - anc + 90, 180) - 90;
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
info = struct('dataset', dataset, 'anchor_axis', anc, 'n_h_simple', nA, 'n_complex_vOS', sum(cls == "complex_vOS"), ...
    'n_orth_simple', sum(cls == "orth_simple"), 'n_other', sum(cls == "other"), 'n_orth_eligible', nOe, ...
    'ambiguous_anchor', amb, 'has_rf', has_rf, 'n_cells', height(T), 'params', o, 'version', 2);
save([od 'types.mat'], 'T', 'info', 'rfimg');
out = info;
if o.do_figs, os_types_figs_claude_auto(T, info, rfimg, Fz.dirs, od, o); end
end
