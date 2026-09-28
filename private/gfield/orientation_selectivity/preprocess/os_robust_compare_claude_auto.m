function P = os_robust_compare_claude_auto(varargin)
% os_robust_compare_claude_auto  Per-cell summaries of os_robust_tuning_claude_auto
% outputs and comparison with GDF curation and the M5 classifier.
% Per dataset, Benjamini-Hochberg FDR over all cell x condition r2 tests (cells with
% >= 1 spike). Per cell:
%   q_min      smallest q over conditions
%   k_sig      conditions with q < alpha and r2 > r1 (orientation > direction
%              component); n_cond conditions; frac_sig = k_sig / n_cond
%   (v2, after QC) r2 and its bootstrap bound are biased upward at low spike counts, so
%   r2x = r2 - r2_null and r2_lbx = r2_lb - r2_null (bias-corrected); conditions with
%   < min_spk (50) spikes are ineligible for 'best' and for k_sig.
%   best cond  = argmax r2_lbx over eligible conditions
%   r2_best, r2_lb_best, r2_lbx_best, r2x_best, r1_best, r1_lb_best, axis_sd_best, rel_best, nspk_best
%   axis_consistency: axial resultant of preferred axes across significant conditions
% Labels: gdf = definite / maybe (os_cells_master_table old_*), gate_rejected
% (passed find_OS_RGCs gate, unlisted), background (neither); new datasets 'new'.
% Output: claude_robust/robust_cells.mat (P) and .csv.
% 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('alpha', 0.05); p.addParameter('min_spk', 50); p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
d = dir([oa 'claude_robust/*_robust.mat']);
Sm = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); OS = Sm.OS;
Q = load([oa 'claude_features/curation_classifier_results.mat'], 'O'); Oc = Q.O;
parts = {};
for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name), 'R'); R = L.R; n = numel(R.cell_ids); nc = size(R.r2, 2);
    pv = R.p_r2(:); ok = ~isnan(pv); qv = nan(size(pv)); qv(ok) = bh(pv(ok)); q = reshape(qv, n, nc);
    sig = q < o.alpha & R.r2 > R.r1;
    ex = R.r2_lb - R.r2_null; ex(R.nspk < o.min_spk) = NaN;   % bias-corrected lower bound; sparse conditions ineligible
    [exb, jb] = max(ex, [], 2); jb(isnan(exb)) = 1; ix = sub2ind([n nc], (1:n)', jb);
    sig = sig & R.nspk >= o.min_spk;
    ac = nan(n, 1);
    for c = 1:n, s = sig(c, :); if any(s), ac(c) = abs(mean(exp(2i * deg2rad(R.pref_axis(c, s))))); end, end
    t = table(repmat(string(R.dataset), n, 1), R.cell_ids(:), repmat(numel(R.dirs), n, 1), repmat(nc, n, 1), ...
        min(q, [], 2), sum(sig, 2), sum(sig, 2) / nc, R.r2(ix), R.r2_lb(ix), exb, R.r2(ix) - R.r2_null(ix), R.r1(ix), R.r1_lb(ix), R.axis_sd(ix), ...
        R.rel(ix), R.nspk(ix), R.pref_axis(ix), ac, R.cond_sp(jb)', R.cond_tp(jb)', ...
        'VariableNames', {'dataset', 'cell_id', 'n_dirs', 'n_cond', 'q_min', 'k_sig', 'frac_sig', 'r2_best', 'r2_lb_best', 'r2_lbx_best', 'r2x_best', ...
        'r1_best', 'r1_lb_best', 'axis_sd_best', 'rel_best', 'nspk_best', 'pref_axis_best', 'axis_consistency', 'best_sp', 'best_tp'});
    parts{end + 1} = t; %#ok<AGROW>
end
P = vertcat(parts{:});
P.label = repmat("background", height(P), 1);
key = P.dataset + "_" + string(P.cell_id);
kO = string(OS.dataset) + "_" + string(OS.cell_id); [tf, lo] = ismember(key, kO);
st = strings(height(P), 1); st(tf) = OS.status(lo(tf));
Ft = load([oa 'claude_features/feature_table.mat'], 'T'); Ft = Ft.T;
kT = string(Ft.dataset) + "_" + string(Ft.cell_id); [tc, lt] = ismember(key, kT);
P.gate_pass = tc; P.osi_best_old = nan(height(P), 1); P.osi_best_old(tc) = Ft.osi_best(lt(tc));
kC = string(Oc.dataset) + "_" + string(Oc.cell_id); [tq, lc] = ismember(key, kC);
P.P_def_lodo = nan(height(P), 1); P.P_def_lodo(tq) = Oc.P_def_lodo(lc(tq));
newds = unique(string(Ft.dataset(Ft.is_new))); isnew = ismember(P.dataset, newds);
P.label(tc & ~isnew) = "gate_rejected";
P.label(st == "old_definite") = "definite"; P.label(st == "old_maybe") = "maybe";
P.label(isnew) = "new"; P.label(isnew & st == "new_definite") = "new_definite"; P.label(isnew & st == "new_maybe") = "new_maybe";
P.Properties.UserData = o;
save([oa 'claude_robust/robust_cells.mat'], 'P'); writetable(P, [oa 'claude_robust/robust_cells.csv']);
end
function q = bh(p)
[ps, i] = sort(p); m = numel(p); qs = ps .* m ./ (1:m)'; qs = flipud(cummin(flipud(qs))); qs(qs > 1) = 1;
q = zeros(size(p)); q(i) = qs;
end
