function P = os_per_trial_pairing_check_claude_auto(dl, varargin)
% os_per_trial_pairing_check_claude_auto  Per-trial check of trial->condition
% pairing, to catch errors confined to a few trials or one stretch of the
% recording (which a dataset-wide median can miss).
%
% For each trial t: the population response (all cells, 0.25 s bins, 0-8 s,
% each cell's trial row z-scored to unit norm) is correlated with
%   (a) the leave-one-out mean of its ASSIGNED condition, and
%   (b) the mean of every other condition sharing all non-DIRECTION params.
% win(t) = (a) > max(b). Correct pairing gives win ~ 1 everywhere; a local
% misalignment appears as a run of losses. Reported: overall win fraction,
% minimum win fraction in any 20-trial window, and that window's start.
%
% 2026-09-27 GDF + Claude
p = inputParser; p.addParameter('win', 20, @isnumeric); p.parse(varargin{:});
datarun = load_data(dl.grating_datapath); datarun = load_neurons(datarun);
datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
S = datarun.stimulus; on = S.triggers(:)'; lab = S.trial_list(:)'; nt = numel(on);
ncell = numel(datarun.spikes); be = 0:0.25:8; nb = numel(be) - 1;
edges = reshape(on + be', 1, []); keepbin = true(1, numel(edges) - 1); keepbin(numel(be):numel(be):end) = false;
Zm = zeros(nt, nb * ncell);
for c = 1:ncell
    h = histcounts(datarun.spikes{c}, edges); X = reshape(h(keepbin), nb, nt)';
    X = X - mean(X, 2); nrm = sqrt(sum(X.^2, 2)); nrm(nrm == 0) = 1;
    Zm(:, (c-1)*nb + (1:nb)) = X ./ nrm;
end
ncomb = numel(S.combinations);
fn = setdiff(fieldnames(S.combinations), {'DIRECTION'});
key = arrayfun(@(q) strjoin(cellfun(@(f) mat2str(S.combinations(q).(f)), fn, 'UniformOutput', false)', '|'), 1:ncomb, 'UniformOutput', false);
[~, ~, grp] = unique(key);
G = sparse(lab, 1:nt, 1, ncomb, nt); Msum = G * Zm; ncnt = full(sum(G, 2));
pc = @(a, b) sum((a - mean(a)) .* (b - mean(b))) / sqrt(sum((a - mean(a)).^2) * sum((b - mean(b)).^2));
win = false(1, nt); ra = nan(1, nt); rb = nan(1, nt);
for t = 1:nt
    q = lab(t); x = Zm(t, :);
    loo = (Msum(q, :) - x) / (ncnt(q) - 1);
    ra(t) = pc(x, loo);
    others = find(grp == grp(q)); others(others == q) = [];
    rb(t) = max(arrayfun(@(o) pc(x, Msum(o, :) / ncnt(o)), others));
    win(t) = ra(t) > rb(t);
end
P.dataset = regexp(dl.grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once');
P.win = win; P.r_assigned = ra; P.r_best_other = rb;
P.win_frac = mean(win);
w = p.Results.win; mv = movmean(double(win), w, 'Endpoints', 'discard');
[P.min_window_frac, P.min_window_start] = min(mv);
P.losing_trials = find(~win);
end
