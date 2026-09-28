function F = os_extract_features_claude_auto(dl, varargin)
% os_extract_features_claude_auto  Per-cell grating-tuning features for one dataset.
%
% usage: F = os_extract_features_claude_auto(dl)   % dl = one registry entry
%
% Loads exactly as the curation batch does (load_data, load_neurons,
% os_trigger_check), then for every cell and every SP x TP condition:
%   osi, dsi, corr   - find_OS_RGCs outputs (the pipeline's threshold metrics)
%   F0, F1n, F2n     - Fourier components of the direction tuning curve
%                      (mean count; |1st| and |2nd| harmonic / (N*F0)); F2n is
%                      an orientation index and F1n a direction index that are
%                      comparable between 8- and 12-direction stimuli
%   pref_ori         - preferred drift-direction AXIS from angle(F2)/2, 0-180 deg
%   pref_dir         - preferred drift direction from angle(F1), 0-360 deg
%   rate             - mean rate (Hz) over directions, 0-8 s window
%   split_half       - corr of tuning curves from odd vs even repetitions
%   tuning           - summed spike counts per direction (cell array per condition)
% Directions are in stimulus coordinates (DIRECTION field of the stim file).
%
% 2026-09-27 GDF + Claude
p = inputParser; p.addParameter('win', 8, @isnumeric); p.parse(varargin{:});
datarun = load_data(dl.grating_datapath);
datarun = load_neurons(datarun);
datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
S = datarun.stimulus;
[osi, dsi, crr] = find_OS_RGCs(datarun, 'all', 'plot_tuning', false, 'verbose', false, 'pause', false);

on = S.triggers(:)'; lab = S.trial_list(:)'; nt = numel(on);
ncell = numel(datarun.spikes);
edges = reshape([on; on + p.Results.win], 1, []);
C = zeros(ncell, nt);
for c = 1:ncell
    h = histcounts(datarun.spikes{c}, edges);
    C(c, :) = h(1:2:end);
end
sps = S.params.SPATIAL_PERIOD; tps = S.params.TEMPORAL_PERIOD; dirs = S.params.DIRECTION;
nd = numel(dirs); th = deg2rad(dirs(:));
comb_dir = arrayfun(@(x) x.DIRECTION, S.combinations);
comb_sp  = arrayfun(@(x) x.SPATIAL_PERIOD, S.combinations);
comb_tp  = arrayfun(@(x) x.TEMPORAL_PERIOD, S.combinations);
z = zeros(ncell, numel(sps), numel(tps));
F = struct('F0', z, 'F1n', z, 'F2n', z, 'pref_ori', z, 'pref_dir', z, 'rate', z, 'split_half', z);
F.tuning = cell(numel(sps), numel(tps));
for a = 1:numel(sps)
    for b = 1:numel(tps)
        T = zeros(ncell, nd); To = T; Te = T;
        for k = 1:nd
            q = find(comb_dir == dirs(k) & comb_sp == sps(a) & comb_tp == tps(b), 1);
            tr = find(lab == q);
            T(:, k)  = sum(C(:, tr), 2);
            To(:, k) = sum(C(:, tr(1:2:end)), 2);
            Te(:, k) = sum(C(:, tr(2:2:end)), 2);
        end
        f0 = mean(T, 2); f1 = T * exp(1i * th); f2 = T * exp(2i * th);
        F.F0(:, a, b)  = f0;
        F.F1n(:, a, b) = abs(f1) ./ (nd * f0);
        F.F2n(:, a, b) = abs(f2) ./ (nd * f0);
        F.pref_ori(:, a, b) = mod(rad2deg(angle(f2)) / 2, 180);
        F.pref_dir(:, a, b) = mod(rad2deg(angle(f1)), 360);
        F.rate(:, a, b) = f0 / S.repetitions / p.Results.win;
        sh = nan(ncell, 1);
        for c = 1:ncell
            if std(To(c, :)) > 0 && std(Te(c, :)) > 0
                r = corrcoef(To(c, :), Te(c, :)); sh(c) = r(1, 2);
            end
        end
        F.split_half(:, a, b) = sh;
        F.tuning{a, b} = T;
    end
end
F.osi = osi; F.dsi = dsi; F.corr = crr;
F.cell_ids = datarun.cell_ids(:)';
F.dataset = regexp(dl.grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once');
F.dirs = dirs; F.sps = sps; F.tps = tps; F.reps = S.repetitions; F.n_trials = nt;
F.grating_datapath = dl.grating_datapath; F.stimulus_path = dl.stimulus_path;
end
