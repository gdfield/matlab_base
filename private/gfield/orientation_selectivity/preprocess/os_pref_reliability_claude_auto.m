function Rl = os_pref_reliability_claude_auto(dataset, cell_ids, best_sp, best_tp, varargin)
% os_pref_reliability_claude_auto  Trial-to-trial reliability of the response in the
% preferred direction, best condition (best_sp, best_tp from os_cells_master_table).
% Preferred direction = direction with the largest mean spike count (1-8 s) in that
% condition. For each repetition, PSTH in bin_s (0.05 s) bins over [t0, t1] (1-8 s,
% excluding the onset transient); reliability = mean pairwise Pearson correlation
% between repetitions (NaN if < 2 repetitions with spikes). Also returns the
% count-based Fano factor at the preferred direction.
% Loading as in os_extract_features_claude_auto (load_data, load_neurons,
% os_trigger_check). 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('bin_s', 0.05); p.addParameter('t0', 1); p.addParameter('t1', 8); p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
reg = load([oa 'claude_features/registry_used.mat']); dl = reg.allds(strcmp(reg.nm, dataset));
dr = load_data(dl.grating_datapath); dr = load_neurons(dr); dr = os_trigger_check(dr, dl.stimulus_path, dl.trigger_interval);
S = dr.stimulus; on = S.triggers(:)'; lab = S.trial_list(:)';
cd_ = arrayfun(@(x) x.DIRECTION, S.combinations); cs = arrayfun(@(x) x.SPATIAL_PERIOD, S.combinations);
ct = arrayfun(@(x) x.TEMPORAL_PERIOD, S.combinations); dirs = S.params.DIRECTION;
e = o.t0:o.bin_s:o.t1; n = numel(cell_ids); z = nan(n, 1);
Rl = table(repmat(string(dataset), n, 1), cell_ids(:), z, z, z, z, 'VariableNames', {'dataset', 'cell_id', 'pref_dir', 'rel_pref', 'fano_pref', 'n_rep'});
for i = 1:n
    ri = find(dr.cell_ids == cell_ids(i)); if isempty(ri), continue; end
    sp = dr.spikes{ri}; mc = nan(numel(dirs), 1); tr = cell(numel(dirs), 1);
    for k = 1:numel(dirs)
        q = find(cd_ == dirs(k) & cs == best_sp(i) & ct == best_tp(i), 1); if isempty(q), continue; end
        tr{k} = find(lab == q);
        mc(k) = mean(arrayfun(@(t) sum(sp >= on(t) + o.t0 & sp < on(t) + o.t1), tr{k}));
    end
    if all(isnan(mc)), continue; end
    [~, kb] = max(mc); t_ = tr{kb}; H = zeros(numel(t_), numel(e) - 1);
    for r = 1:numel(t_), H(r, :) = histcounts(sp - on(t_(r)), e); end
    cnt = sum(H, 2); Rl.pref_dir(i) = dirs(kb); Rl.n_rep(i) = numel(t_);
    if mean(cnt) > 0, Rl.fano_pref(i) = var(cnt) / mean(cnt); end
    ok = std(H, 0, 2) > 0;
    if sum(ok) >= 2, c = corrcoef(H(ok, :)'); Rl.rel_pref(i) = mean(c(triu(true(sum(ok)), 1))); end
end
end
