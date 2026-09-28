function R = os_robust_tuning_claude_auto(dl, varargin)
% os_robust_tuning_claude_auto  Trial-level orientation-tuning statistics for every
% cell and every SP x TP condition of one dataset (GDF 2026-09-28: statistical test,
% bootstrap lower bounds, condition specificity, direction-sampling-independent metrics).
%
% Loading and trial parsing exactly as os_extract_features_claude_auto (load_data,
% load_neurons, os_trigger_check); spike counts in [onset, onset + win) (win 8 s).
% Per cell and condition, with m_k = mean count at direction theta_k:
%   r2 = |sum_k m_k exp(2i theta_k)| / sum_k m_k   (vector-sum orientation index;
%        = F2n; 1 - r2 = circular variance of the axial tuning; for evenly spaced
%        directions it does not depend on the number of directions, 8 vs 12)
%   r1 = same with exp(1i theta_k) (direction index); pref_axis = angle(F2)/2.
%   p_r2, p_r1: permutation test, n_perm (1000) random reassignments of the
%        condition's trials to directions (preserving the number of trials per
%        direction); p = (1 + #null >= observed) / (1 + n_perm).
%   r2_lb, r1_lb: 5th percentile of n_boot (1000) bootstrap replicates (repetitions
%        resampled with replacement within each direction), so low spike counts and
%        trial-to-trial variability lower the bound.
%   axis_sd: circular SD (deg, axial) of the bootstrap preferred axis.
%   r2_null, r1_null: mean of the permutation null (the index expected from noise
%        alone at this cell's spike count; r2 and its bootstrap bound are biased
%        upward by about this amount at low counts; QC 2026-09-28).
%   rel: mean pairwise correlation across repetitions of the direction tuning vector.
%   nspk: total spikes in the condition.
% 'keep_ids': cell ids whose permutation-null and bootstrap distributions (and the
% resampling counts and trial counts) are returned in R.D, for illustration.
% Saves os_analysis/claude_robust/<dataset>_robust.mat (struct R) unless 'save' is false.
% 'calib_shuffle' (calibration only): trials are shuffled across directions before
% analysis, so p values should be uniform.
% 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('win', 8); p.addParameter('n_perm', 1000); p.addParameter('n_boot', 1000);
p.addParameter('seed', 1); p.addParameter('calib_shuffle', false); p.addParameter('save', true); p.addParameter('keep_ids', []);
p.parse(varargin{:}); o = p.Results;
rng(o.seed, 'twister'); D = struct();
datarun = load_data(dl.grating_datapath); datarun = load_neurons(datarun);
datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
S = datarun.stimulus; on = S.triggers(:)'; lab = S.trial_list(:)'; nt = numel(on);
ncell = numel(datarun.spikes); edges = reshape([on; on + o.win], 1, []);
if any(diff(edges) <= 0), error('trial windows overlap (onsets closer than win)'); end
C = zeros(ncell, nt);
for c = 1:ncell, h = histcounts(datarun.spikes{c}, edges); C(c, :) = h(1:2:end); end
sps = S.params.SPATIAL_PERIOD; tps = S.params.TEMPORAL_PERIOD; dirs = S.params.DIRECTION; nd = numel(dirs);
cd_ = arrayfun(@(x) x.DIRECTION, S.combinations); cs = arrayfun(@(x) x.SPATIAL_PERIOD, S.combinations);
ct = arrayfun(@(x) x.TEMPORAL_PERIOD, S.combinations);
nc = numel(sps) * numel(tps); z = nan(ncell, nc);
R = struct('r2', z, 'r1', z, 'p_r2', z, 'p_r1', z, 'r2_lb', z, 'r1_lb', z, 'pref_axis', z, 'axis_sd', z, ...
    'rel', z, 'nspk', z, 'rate', z, 'r2_null', z, 'r1_null', z);
cond_sp = nan(1, nc); cond_tp = cond_sp; j = 0;
for a = 1:numel(sps)
    for b = 1:numel(tps)
        j = j + 1; cond_sp(j) = sps(a); cond_tp(j) = tps(b);
        tr = cell(1, nd); dk = [];
        for k = 1:nd
            q = find(cd_ == dirs(k) & cs == sps(a) & ct == tps(b), 1); tr{k} = find(lab == q); dk = [dk, k * ones(1, numel(tr{k}))]; %#ok<AGROW>
        end
        ti = [tr{:}]; Cc = C(:, ti); nrep = cellfun(@numel, tr);
        if o.calib_shuffle, Cc = Cc(:, randperm(numel(ti))); end   % calibration: destroy direction information
        th = deg2rad(dirs(dk)); th = th(:);
        % observed (mean per direction -> weight 1/nrep)
        w = 1 ./ nrep(dk); w = w(:);
        e2 = w .* exp(2i * th); e1 = w .* exp(1i * th); s0 = Cc * w;
        f2 = Cc * e2; f1 = Cc * e1;
        R.r2(:, j) = abs(f2) ./ s0; R.r1(:, j) = abs(f1) ./ s0; R.pref_axis(:, j) = mod(rad2deg(angle(f2)) / 2, 180);
        R.nspk(:, j) = sum(Cc, 2); R.rate(:, j) = s0 / nd / o.win;
        % permutation null: shuffle direction labels over the condition's trials
        P2 = zeros(numel(ti), o.n_perm); P1 = P2; P0 = P2;
        for i = 1:o.n_perm
            pk = dk(randperm(numel(dk))); tp_ = deg2rad(dirs(pk)); wp = 1 ./ nrep(pk);
            P0(:, i) = wp(:); P2(:, i) = wp(:) .* exp(2i * tp_(:)); P1(:, i) = wp(:) .* exp(1i * tp_(:));
        end
        sp0 = Cc * P0; n2 = abs(Cc * P2) ./ sp0; n1 = abs(Cc * P1) ./ sp0;
        tol = 1e-12;   % floating-point ties (|exp(2i theta)| ~= 1 - eps)
        R.p_r2(:, j) = (1 + sum(n2 >= R.r2(:, j) - tol, 2)) / (1 + o.n_perm);
        if ~isempty(o.keep_ids), kc = find(ismember(datarun.cell_ids(:), o.keep_ids(:))); D.null_r2{j} = n2(kc, :); D.Cc{j} = Cc(kc, :); D.dk{j} = dk; D.kc_ids = datarun.cell_ids(kc); end
        R.p_r1(:, j) = (1 + sum(n1 >= R.r1(:, j) - tol, 2)) / (1 + o.n_perm);
        R.r2_null(:, j) = mean(n2, 2); R.r1_null(:, j) = mean(n1, 2);
        z0 = s0 <= 0; R.p_r2(z0, j) = NaN; R.p_r1(z0, j) = NaN;
        % bootstrap over repetitions within direction
        B2 = zeros(numel(ti), o.n_boot); B0 = B2; B1 = B2; off = [0 cumsum(nrep)];
        for i = 1:o.n_boot
            cnt = zeros(numel(ti), 1);
            for k = 1:nd, r_ = randi(nrep(k), nrep(k), 1); cnt(off(k) + (1:nrep(k))) = accumarray(r_, 1, [nrep(k) 1]); end
            wb = cnt .* w; B0(:, i) = wb; B2(:, i) = wb .* exp(2i * th); B1(:, i) = wb .* exp(1i * th);
        end
        sb = Cc * B0; b2 = Cc * B2; b1 = Cc * B1;
        if ~isempty(o.keep_ids), D.boot_r2{j} = abs(b2(kc, :)) ./ sb(kc, :); D.boot_cnt{j} = B0 ./ w; D.nrep{j} = nrep; end
        R.r2_lb(:, j) = prctile(abs(b2) ./ sb, 5, 2); R.r1_lb(:, j) = prctile(abs(b1) ./ sb, 5, 2);
        ang = angle(b2); R.axis_sd(:, j) = rad2deg(sqrt(-2 * log(abs(mean(exp(1i * ang), 2))))) / 2;
        % reliability across repetitions (tuning vector per repetition index)
        nr = min(nrep); V = zeros(ncell, nd, nr);
        for k = 1:nd, V(:, k, :) = reshape(C(:, tr{k}(1:nr)), ncell, 1, nr); end
        rl = nan(ncell, 1);
        for c = 1:ncell
            X = squeeze(V(c, :, :)); if nr < 2 || any(std(X) == 0), continue; end
            cc = corrcoef(X); rl(c) = mean(cc(triu(true(nr), 1)));
        end
        R.rel(:, j) = rl;
    end
end
R.cell_ids = datarun.cell_ids(:)'; R.cond_sp = cond_sp; R.cond_tp = cond_tp; R.dirs = dirs;
R.reps = S.repetitions; R.dataset = regexp(dl.grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once');
R.params = o; if ~isempty(o.keep_ids), R.D = D; end
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/claude_robust/';
if o.save, if ~exist(oa, 'dir'), mkdir(oa); end, save([oa R.dataset '_robust.mat'], 'R'); end
end
