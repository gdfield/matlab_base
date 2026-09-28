function A = ei_axon_direction_claude_auto(ei, pos, varargin)
% ei_axon_direction_claude_auto  Axon propagation direction from one electrical image.
%
% usage: A = ei_axon_direction_claude_auto(ei, pos)
%   ei  : n_electrodes x n_samples (Vision EI, 20 kHz)
%   pos : n_electrodes x 2 electrode positions (array coordinates, um)
% Soma electrode = largest peak-to-peak amplitude. Axon electrodes = electrodes
% whose negative peak occurs >= min_lag samples after the soma's, that lie
% >= min_dist um away, and whose amplitude exceeds max(frac*soma amp, abs_thr).
% Direction = amplitude-weighted mean of the unit vectors soma -> axon electrode.
% Returns angle (deg, 0-360, array coords; direction of propagation, i.e. toward
% the optic disc), resultant length R (0-1, spread of axon electrodes), n_axon,
% speed (um/ms, slope of distance vs lag; positive for propagation away from soma).
% 2026-09-27 GDF + Claude
p = inputParser;
p.addParameter('frac', 0.05); p.addParameter('abs_thr', 1); p.addParameter('min_lag', 2);
p.addParameter('min_dist', 60); p.addParameter('fs', 20000); p.parse(varargin{:}); o = p.Results;
A = struct('angle', NaN, 'R', NaN, 'n_axon', 0, 'speed', NaN, 'soma_xy', [NaN NaN]);
if isempty(ei), return; end
amp = max(ei, [], 2) - min(ei, [], 2);
[~, tt] = min(ei, [], 2);
[as, s] = max(amp); A.soma_xy = pos(s, :);
d = pos - pos(s, :); dist = sqrt(sum(d.^2, 2));
ax = tt >= tt(s) + o.min_lag & dist >= o.min_dist & amp >= max(o.frac * as, o.abs_thr);
A.n_axon = sum(ax);
if A.n_axon < 3, return; end
u = d(ax, :) ./ dist(ax); w = amp(ax);
v = sum(u .* w, 1) / sum(w);
A.R = norm(v); A.angle = mod(atan2d(v(2), v(1)), 360);
lag_ms = (tt(ax) - tt(s)) / o.fs * 1000;
if numel(unique(lag_ms)) > 1, b = polyfit(lag_ms, dist(ax), 1); A.speed = b(1); end
end
