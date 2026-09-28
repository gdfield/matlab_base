function R = os_pairing_shift_test_claude_auto(dl, varargin)
% os_pairing_shift_test_claude_auto  Test that grating trials are paired with
% the correct stimulus conditions, for one dataset registry entry.
%
% usage: R = os_pairing_shift_test_claude_auto(dl)
%   dl: one entry of new_os_datasets() / os_datasets() (grating_datapath,
%       stimulus_path, trigger_interval). Loads through os_trigger_check,
%       exactly as the curation batch does.
%
% Statistic: for each cell, within-condition reliability = mean pairwise
% Pearson correlation between the binned spike trains (0.25 s bins, 0-8 s
% after trial onset) of repeats of the same stimulus condition, averaged over
% conditions with >= 3 non-silent repeats and >= 10 spikes. The dataset
% statistic is the median over cells.
%
% Tests:
%   1) shift profile: re-pair trial onset t with file trial t+k, k = -5..5.
%      Correct pairing should give the maximum at k = 0.
%   2) permutation test: shuffle DIRECTION labels among trials that share all
%      other parameters (SP, TP, RGB, ...), n_perm times (default 1000).
%      p = (1 + #null >= observed) / (n_perm + 1).
%
% Mean pairwise correlation is computed exactly from z-scored trial vectors:
% for unit-norm, zero-mean rows z_i in a group of n,
%   mean_{i<j} corr(z_i,z_j) = (||sum_i z_i||^2 - n) / (n(n-1)).
% Silent (zero-variance) trials are excluded, matching corrcoef NaN handling.
%
% 2026-09-27 GDF + Claude

p = inputParser;
p.addParameter('n_perm', 1000, @isnumeric);
p.addParameter('shifts', -5:5, @isnumeric);
p.addParameter('bin', 0.25, @isnumeric);
p.addParameter('win', 8, @isnumeric);
p.addParameter('seed', 1, @isnumeric);
p.parse(varargin{:});
rng(p.Results.seed);

datarun = load_data(dl.grating_datapath);
datarun = load_neurons(datarun);
datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
S = datarun.stimulus;
on = S.triggers(:)';
lab = S.trial_list(:)';
nt = numel(on);
ncell = numel(datarun.spikes);

% --- binned responses: Z (nt x nb x ncell) z-scored rows, V valid, C counts
be = 0:p.Results.bin:p.Results.win; nb = numel(be) - 1;
edges = reshape(on + be', 1, []);           % per-trial edge blocks
keepbin = true(1, numel(edges) - 1);
keepbin(numel(be):numel(be):end) = false;   % drop gap bins between trials
Z = zeros(nt, nb, ncell); V = false(nt, ncell); C = zeros(nt, ncell);
for c = 1:ncell
    h = histcounts(datarun.spikes{c}, edges);
    X = reshape(h(keepbin), nb, nt)';
    C(:, c) = sum(X, 2);
    Xc = X - mean(X, 2);
    nrm = sqrt(sum(Xc.^2, 2));
    V(:, c) = nrm > 0;
    Xc(V(:, c), :) = Xc(V(:, c), :) ./ nrm(V(:, c));
    Xc(~V(:, c), :) = 0;
    Z(:, :, c) = Xc;
end
Zm = reshape(Z, nt, nb * ncell);

% --- condition grouping for the direction shuffle
ncomb = numel(S.combinations);
fn = setdiff(fieldnames(S.combinations), {'DIRECTION'});
key = cell(ncomb, 1);
for q = 1:ncomb
    parts = cellfun(@(f) mat2str(S.combinations(q).(f)), fn, 'UniformOutput', false);
    key{q} = strjoin(parts', '|');
end
[~, ~, grp_of_comb] = unique(key);

stat = @(labv, rows) score(labv, rows);

    function m = score(labv, rows)
        G = sparse(labv, 1:numel(rows), 1, ncomb, numel(rows));
        Sm = G * Zm(rows, :);
        sumsq = squeeze(sum(reshape(Sm, ncomb, nb, ncell).^2, 2));
        n = G * double(V(rows, :));
        sp = G * C(rows, :);
        r = (sumsq - n) ./ (n .* (n - 1));
        ok = n >= 3 & sp >= 10;
        r(~ok) = NaN;
        cellr = mean(r, 1, 'omitnan');
        m = median(cellr(~isnan(cellr)));
    end

% --- observed and shift profile
R.dataset = regexp(dl.grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once');
R.n_trials = nt; R.n_cells = ncell; R.n_reps = S.repetitions; R.n_combos = ncomb;
R.shifts = p.Results.shifts;
R.shift_profile = nan(1, numel(R.shifts));
for si = 1:numel(R.shifts)
    k = R.shifts(si);
    t = 1:nt; t = t(t + k >= 1 & t + k <= nt);
    R.shift_profile(si) = stat(lab(t + k), t);
end
R.observed = R.shift_profile(R.shifts == 0);
[~, im] = max(R.shift_profile); R.best_shift = R.shifts(im);
R.max_other_shift = max(R.shift_profile(R.shifts ~= 0));

% --- permutation null: shuffle DIRECTION within groups sharing all other params
tg = grp_of_comb(lab)';
% map (group, direction) -> combination index
dirv = arrayfun(@(x) x.DIRECTION, S.combinations);
null = nan(1, p.Results.n_perm);
for pi = 1:p.Results.n_perm
    newlab = lab;
    for g = unique(tg)
        idx = find(tg == g);
        newlab(idx) = lab(idx(randperm(numel(idx))));
    end
    null(pi) = stat(newlab, 1:nt);
end
R.null_mean = mean(null); R.null_p95 = prctile(null, 95); R.null_max = max(null);
R.p_perm = (1 + sum(null >= R.observed)) / (p.Results.n_perm + 1);
R.n_dir_per_group = numel(unique(dirv));
end
