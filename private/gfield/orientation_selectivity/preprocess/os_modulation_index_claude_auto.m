function MI = os_modulation_index_claude_auto(dl, varargin)
% os_modulation_index_claude_auto  Simple/complex (linear/nonlinear) index per cell.
%
% For every cell: PSTH (10 ms bins, 1-8 s after trial onset, skipping the onset
% transient) for each SP x TP condition and direction. Modulation index
% MI = 2|F1| / F0 of the PSTH at the grating temporal frequency (TF), reported at
% the cell's best condition (argmax F2n * max(corr,0), as in the feature table)
% and preferred direction (max spike count). Simple (hOS) cells follow the
% grating cycle (high MI); complex (vOS) cells do not (low MI).
% TF = frame_rate / TP, with the display frame rate chosen from the data (see below;
% 120 Hz in 2012 recordings, 60 Hz later), rather than assumed. This avoids a wrong
% TF silently driving MI toward zero.
% 2026-09-27 GDF + Claude
p = inputParser; p.addParameter('bin', 0.01); p.addParameter('t0', 1); p.addParameter('t1', 8);
p.parse(varargin{:}); o = p.Results;
datarun = load_data(dl.grating_datapath); datarun = load_neurons(datarun);
datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
S = datarun.stimulus; on = S.triggers(:)'; lab = S.trial_list(:)'; nt = numel(on);
be = o.t0:o.bin:o.t1; nb = numel(be) - 1; tc = be(1:end-1) + o.bin/2;
edges = reshape(on + be', 1, []); keep = true(1, numel(edges) - 1); keep(numel(be):numel(be):end) = false;
sps = S.params.SPATIAL_PERIOD; tps = S.params.TEMPORAL_PERIOD; dirs = S.params.DIRECTION; nd = numel(dirs);
cd = arrayfun(@(x) x.DIRECTION, S.combinations); cs = arrayfun(@(x) x.SPATIAL_PERIOD, S.combinations);
ct = arrayfun(@(x) x.TEMPORAL_PERIOD, S.combinations);
nc = numel(datarun.spikes); ncond = numel(sps) * numel(tps);
% trial index sets per (condition k, direction j); k runs SP fastest
trials = cell(ncond, nd); tpk = zeros(1, ncond); k = 0;
for b = 1:numel(tps)
    for a = 1:numel(sps)
        k = k + 1; tpk(k) = b;
        for j = 1:nd
            trials{k, j} = find(lab == find(cd == dirs(j) & cs == sps(a) & ct == tps(b), 1));
        end
    end
end
Pc = zeros(nc, ncond, nd, nb, 'single');
pop = zeros(numel(tps), nb);
for c = 1:nc
    h = histcounts(datarun.spikes{c}, edges); X = reshape(h(keep), nb, nt)';
    for k = 1:ncond
        for j = 1:nd
            v = sum(X(trials{k, j}, :), 1); Pc(c, k, j, :) = v; pop(tpk(k), :) = pop(tpk(k), :) + v;
        end
    end
end
% Display frame rate chosen among plausible values from the responses themselves:
% for each cell (>= 200 spikes) and TP, the PSTH summed over that TP's trials;
% spectral SNR at TF = rate/TP = amplitude at TF / median amplitude within
% +/-0.5 Hz (excluding +/-0.1 Hz). Score = median SNR over cells, averaged over
% TPs; one rate per dataset. (Checked: 480- and 540-frame stimulus files -> 60 Hz,
% 960-frame files -> 120 Hz.)
fr_cand = [30 60 75 85 120]; fgrid = 0.05:0.02:12; E = exp(-2i * pi * fgrid(:) * tc);
amp = nan(numel(fr_cand), numel(tps));
for b = 1:numel(tps)
    Y = squeeze(sum(sum(Pc(:, tpk == b, :, :), 2), 3));   % nc x nb
    if nc == 1, Y = Y(:)'; end
    ok = sum(Y, 2) >= 200; vals = nan(nc, numel(fr_cand));
    A = abs(E * double(Y(ok, :) - mean(Y(ok, :), 2))');   % nf x n_ok
    for m = 1:numel(fr_cand)
        f0 = fr_cand(m) / tps(b); if f0 > 11.5 || f0 < 0.3, continue; end
        [~, i0] = min(abs(fgrid - f0)); nbh = abs(fgrid - f0) <= 0.5 & abs(fgrid - f0) > 0.1;
        vals(ok, m) = (A(i0, :) ./ median(A(nbh, :), 1))';
    end
    amp(:, b) = median(vals, 1, 'omitnan')';
end
[~, im] = max(mean(amp, 2, 'omitnan')); frame_rate = fr_cand(im); tf_est = frame_rate ./ tps;
Nd = sum(Pc, 4);                              % nc x ncond x nd spike counts
th = reshape(deg2rad(dirs), 1, 1, nd);
f2 = abs(sum(Nd .* exp(2i * th), 3)) ./ (nd * mean(Nd, 3));
[~, ~, crr3] = find_OS_RGCs(datarun, 'all', 'plot_tuning', false, 'verbose', false, 'pause', false);
crr = reshape(crr3, nc, []);
[~, best] = max(f2 .* max(crr, 0), [], 2);
mi = nan(nc, ncond);
for k = 1:ncond
    w = exp(-2i * pi * tf_est(tpk(k)) * tc);
    [~, jp] = max(squeeze(Nd(:, k, :)), [], 2);
    for c = 1:nc
        ps = double(squeeze(Pc(c, k, jp(c), :)))';
        mi(c, k) = 2 * abs(sum(ps .* w)) / max(sum(ps), eps);
    end
end
MI = table(datarun.cell_ids(:), mi(sub2ind(size(mi), (1:nc)', best)), best, 'VariableNames', {'cell_id', 'mod_index', 'best_cond'});
MI.Properties.UserData = struct('mi_all', mi, 'sps', sps, 'tps', tps, 'tf_est', tf_est, 'frame_rate_est', frame_rate, 'fr_cand', fr_cand, 'fr_amp', amp, ...
    'dataset', regexp(dl.grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once'));
end
