function plot_chirp_psth(datarun, cell_id, varargin)
% PLOT_CHIRP_PSTH   Plot one cell's trial-averaged chirp PSTH with phase regions shaded
%
% usage: plot_chirp_psth(datarun, cell_id, varargin)
%
% arguments:  datarun  - datarun struct with datarun.chirp populated (see LOAD_CHIRP)
%             cell_id  - single cell id to plot
%             varargin - struct or list of optional parameters (see below)
%
% optional params, their default values, and what they specify:
%
% figure    0     figure or axes to plot in, see SET_UP_FIG_OR_AXES
%                     (0 = new figure, -1 = current axes)
%
% See also: LOAD_CHIRP, PLOT_CHIRP_RASTER, PLOT_STA
%
% 2026-07 mcarleton
%

p = inputParser;
p.addParamValue('figure', 0);
p.parse(varargin{:});
params = p.Results;

if ~isfield(datarun,'chirp')
    error('plot_chirp_psth: datarun.chirp not found - run load_chirp(datarun, chirp_path) first.');
end

cell_num = get_cell_indices(datarun, cell_id);
if ~datarun.chirp.matched(cell_num)
    error('plot_chirp_psth: cell %d has no chirp data (no match in the loaded chirp file).', cell_id);
end

psth = datarun.chirp.psth_mean{cell_num};
bin_size = datarun.chirp.params.bin_size_ms;
t_ms = (0:length(psth)-1) * bin_size;

plot_axes = set_up_fig_or_axes(params.figure);
axes(plot_axes);
hold on;

% shade each stimulus phase using the bin ranges loaded by load_chirp
phase_colors = struct( ...
    'step_on',    [0.7 0.7 0.7], ...
    'step_off',   [0.4 0.4 0.4], ...
    'freq_sweep', [0.6 0.7 1.0], ...
    'contrast',   [0.6 1.0 0.6]);

y_max = max(psth) * 1.05 + eps;
phase_names = fieldnames(datarun.chirp.phase);
legend_handles = [];
legend_labels = {};

for pn = 1:length(phase_names)
    name = phase_names{pn};
    se = datarun.chirp.phase.(name);
    t_start = (se(1) - 1) * bin_size;
    t_end = se(2) * bin_size;
    hp = patch([t_start t_end t_end t_start], [0 0 y_max y_max], phase_colors.(name), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4);
    legend_handles(end+1) = hp; %#ok<AGROW>
    legend_labels{end+1} = strrep(name, '_', ' '); %#ok<AGROW>
end

hpsth = plot(t_ms, psth, 'k-', 'LineWidth', 1.5);
legend_handles(end+1) = hpsth;
legend_labels{end+1} = 'PSTH';

xlabel('time from stim onset (ms)');
ylabel('firing rate (Hz)');
title(sprintf('cell %d chirp PSTH (QI = %.2f)', cell_id, datarun.chirp.quality_index(cell_num)));
legend(legend_handles, legend_labels, 'Location', 'northeastoutside');
xlim([0 t_ms(end)]);
ylim([0 y_max]);
hold off;
