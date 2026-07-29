function plot_chirp_raster(datarun, cell_id, varargin)
% PLOT_CHIRP_RASTER   Plot per-trial binned spike counts for one cell as a trial x time heatmap
%
% NOTE: the chirp analysis pipeline only saves per-trial BINNED spike counts
% (spikes_binned), not individual spike timestamps - so this is not a classic
% spike-time dot raster. It plots trial-by-trial binned counts as a heatmap
% (rows = trials, columns = time bins), which is the closest thing available
% from this file. For a true dot raster you'd need the raw datarun.spikes plus
% the chirp stimulus's trigger times, which this file does not provide -
% see TESTRASTER / BROWSE_RASTERS for that approach on other stimulus types.
%
% usage: plot_chirp_raster(datarun, cell_id, varargin)
%
% arguments:  datarun  - datarun struct with datarun.chirp populated (see LOAD_CHIRP)
%             cell_id  - single cell id to plot
%             varargin - struct or list of optional parameters (see below)
%
% optional params, their default values, and what they specify:
%
% figure     0        figure or axes to plot in, see SET_UP_FIG_OR_AXES
% colormap   'hot'     matlab colormap name
%
% See also: LOAD_CHIRP, PLOT_CHIRP_PSTH, TESTRASTER, BROWSE_RASTERS
%
% 2026-07 mcarleton
%

p = inputParser;
p.addParamValue('figure', 0);
p.addParamValue('colormap', 'hot');
p.parse(varargin{:});
params = p.Results;

if ~isfield(datarun,'chirp')
    error('plot_chirp_raster: datarun.chirp not found - run load_chirp(datarun, chirp_path) first.');
end

cell_num = get_cell_indices(datarun, cell_id);
if ~datarun.chirp.matched(cell_num)
    error('plot_chirp_raster: cell %d has no chirp data (no match in the loaded chirp file).', cell_id);
end

trials = datarun.chirp.spikes_binned{cell_num};   % (n_trials, n_bins) raw counts
bin_size = datarun.chirp.params.bin_size_ms;
t_ms = (0:size(trials,2)-1) * bin_size;

plot_axes = set_up_fig_or_axes(params.figure);
axes(plot_axes);

imagesc(t_ms, 1:size(trials,1), trials);
colormap(plot_axes, params.colormap);
colorbar;
xlabel('time from stim onset (ms)');
ylabel('trial #');
title(sprintf('cell %d - per-trial binned spike counts (QI = %.2f)', ...
    cell_id, datarun.chirp.quality_index(cell_num)));
