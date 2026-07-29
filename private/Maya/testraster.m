function [trial_spike_times, firing_rate, bin_centers] = testraster(datarun1, cell_id, varargin)
%TESTRASTER  Plot a trial raster + PSTH for one cell, at one stimulus
%condition or across all conditions.
%
%   testraster(datarun, cell_id)
%   testraster(datarun, cell_id, 'direction', 90)
%   testraster(datarun, cell_id, 'temporal_period', 256)
%   testraster(datarun, cell_id, 'direction', 90, 'temporal_period', 256)
%   testraster(datarun, cell_id, 'direction', 90, 'bin_width', 0.02)
%
%   [trial_spike_times, firing_rate, bin_centers] = testraster(...)
%
% REQUIRED INPUTS
%   datarun   - the datarun struct (after load_data/load_neurons/
%               load_params/load_stim)
%   cell_id   - the RGC cell ID to plot (e.g. datarun.cell_ids(1))
%
% OPTIONAL NAME-VALUE INPUTS
%   'direction'         - stimulus DIRECTION to plot (e.g. 90). If
%                          omitted, all directions are included.
%   'temporal_period'   - stimulus TEMPORAL_PERIOD to plot (e.g. 256).
%                          If omitted, all temporal periods are included.
%   'bin_width'         - PSTH bin width in seconds (default 0.05)
%   'plot'              - true/false, whether to generate a figure
%                          (default true)
%
% BEHAVIOR
%   - If both 'direction' and 'temporal_period' are given, exactly one
%     stimulus combination is matched and a single raster+PSTH figure is
%     plotted.
%   - If either (or both) is omitted, every matching combination is
%     included, and results are shown as a grid of raster subpanels (one
%     per combination) in a single figure.
%
% OUTPUTS
%   For a single matched condition:
%       trial_spike_times - cell array, one entry per trial, of spike
%                            times (s) relative to stimulus onset
%       firing_rate        - PSTH firing rate (spikes/s), pooled across trials
%       bin_centers         - time (s) at the center of each PSTH bin
%   For multiple matched conditions:
%       trial_spike_times - cell array indexed by matched condition,
%                            each entry itself a cell array of per-trial
%                            spike times (as above)
%       firing_rate         - cell array indexed by matched condition
%       bin_centers         - time (s) at the center of each PSTH bin
%                            (same for every condition)

%% parse inputs
p = inputParser;
addParameter(p, 'direction', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'temporal_period', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'bin_width', 0.05, @isnumeric);
addParameter(p, 'plot', true, @islogical);
parse(p, varargin{:});

direction       = p.Results.direction;
temporal_period = p.Results.temporal_period;
bin_width       = p.Results.bin_width;
do_plot         = p.Results.plot;

%% pull needed values out of datarun
stim_dur = 8;   % change here (or promote to an input) if your stimulus duration differs
if isfield(datarun1.stimulus, 'repetitions') && ~isempty(datarun1.stimulus.repetitions)
    num_reps = datarun1.stimulus.repetitions;
else
    % fall back to computing it directly from the trial list, in case
    % this datarun wasn't loaded via get_all_condition_spike_times or
    % repetitions wasn't set for some other reason
    trial_counts = histcounts(datarun1.stimulus.trial_list, 0.5:1:(length(datarun1.stimulus.combinations) + 0.5));
    num_reps = mode(trial_counts);
end

cell_idx = get_cell_indices(datarun1, cell_id);
spike_times = datarun1.spikes{cell_idx};   % all spike times (s) for this cell

%% find stimulus combination(s) matching the requested direction/TP
all_dirs = [datarun1.stimulus.combinations.DIRECTION];
all_tps  = [datarun1.stimulus.combinations.TEMPORAL_PERIOD];

match = true(size(all_dirs));
if ~isempty(direction)
    match = match & (all_dirs == direction);
end
if ~isempty(temporal_period)
    match = match & (all_tps == temporal_period);
end

stim_indices = find(match);

if isempty(stim_indices)
    error('testraster:noMatch', ...
        'No stimulus combination found matching direction=%s, temporal_period=%s.', ...
        mat2str(direction), mat2str(temporal_period));
end

%% PSTH bin edges (shared across conditions)
edges = 0:bin_width:stim_dur;
bin_centers = edges(1:end-1) + bin_width/2;
tick_halfheight = 0.4;

%% -------- CASE 1: exactly one condition matched --------
if isscalar(stim_indices)
    stim_index = stim_indices;
    this_combo = datarun1.stimulus.combinations(stim_index);

    trial_inds = find(datarun1.stimulus.trial_list == stim_index);
    if length(trial_inds) ~= num_reps
        warning('testraster:repCount', ...
            'Found %d trials for this condition, expected num_reps = %d', ...
            length(trial_inds), num_reps);
    end
    trial_start_times = datarun1.stimulus.triggers(trial_inds);
    n_trials = length(trial_inds);

    trial_spike_times = cell(n_trials, 1);
    for rep = 1:n_trials
        t0 = trial_start_times(rep);
        t1 = t0 + stim_dur;
        trial_spike_times{rep} = spike_times(spike_times >= t0 & spike_times < t1) - t0;
    end

    all_spikes = vertcat(trial_spike_times{:});
    counts = histcounts(all_spikes, edges);
    firing_rate = counts / (n_trials * bin_width);

    if do_plot
        figure('Name', sprintf('Cell %d - dir=%g, TP=%g', ...
            cell_id, this_combo.DIRECTION, this_combo.TEMPORAL_PERIOD));

        subplot(2,1,1); hold on;
        for rep = 1:n_trials
            spk = trial_spike_times{rep};
            for s = 1:length(spk)
                line([spk(s) spk(s)], [rep - tick_halfheight, rep + tick_halfheight], ...
                    'Color', 'k', 'LineWidth', 1);
            end
        end
        xlim([0 stim_dur]);
        ylim([0.5 n_trials + 0.5]);
        set(gca, 'YTick', 1:n_trials, 'XTickLabel', []);
        ylabel('Trial');
        title(sprintf('Cell %d, dir=%g, TP=%g (%d repeats)', ...
            cell_id, this_combo.DIRECTION, this_combo.TEMPORAL_PERIOD, n_trials));
        box on;

        subplot(2,1,2);
        bar(bin_centers, firing_rate, 'FaceColor', 'k', 'EdgeColor', 'none', 'BarWidth', 1);
        xlim([0 stim_dur]);
        xlabel('Time from stimulus onset (s)');
        ylabel('Firing rate (sp/s)');
        box on;
    end

    return
end

%% -------- CASE 2: multiple conditions matched -> grid of subpanels --------
n_cond = length(stim_indices);

% order matched conditions by direction (ties keep original relative order)
matched_dirs = all_dirs(stim_indices);
[~, order] = sort(matched_dirs);
stim_indices = stim_indices(order);

trial_spike_times = cell(n_cond, 1);
firing_rate = cell(n_cond, 1);

if do_plot
    n_cols = ceil(sqrt(n_cond));
    n_rows = ceil(n_cond / n_cols);
    figure('Name', sprintf('Cell %d - matched stimulus conditions', cell_id));
    t = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
end

for ii = 1:n_cond
    stim_index = stim_indices(ii);
    this_combo = datarun1.stimulus.combinations(stim_index);

    trial_inds = find(datarun1.stimulus.trial_list == stim_index);
    trial_start_times = datarun1.stimulus.triggers(trial_inds);
    n_trials = length(trial_inds);

    cond_spikes = cell(n_trials, 1);
    for rep = 1:n_trials
        t0 = trial_start_times(rep);
        t1 = t0 + stim_dur;
        cond_spikes{rep} = spike_times(spike_times >= t0 & spike_times < t1) - t0;
    end
    trial_spike_times{ii} = cond_spikes;

    all_spikes = vertcat(cond_spikes{:});
    counts = histcounts(all_spikes, edges);
    firing_rate{ii} = counts / (max(n_trials,1) * bin_width);

    if do_plot
        nexttile; hold on;
        for rep = 1:n_trials
            spk = cond_spikes{rep};
            for s = 1:length(spk)
                line([spk(s) spk(s)], [rep - tick_halfheight, rep + tick_halfheight], ...
                    'Color', 'k', 'LineWidth', 0.75);
            end
        end
        xlim([0 stim_dur]);
        ylim([0.5 num_reps + 0.5]);
        set(gca, 'YTick', [], 'XTick', 0:2:stim_dur, 'FontSize', 7, 'Box', 'on');
        title(sprintf('dir=%g, TP=%g', this_combo.DIRECTION, this_combo.TEMPORAL_PERIOD), ...
            'FontSize', 8, 'FontWeight', 'normal');
    end
end

if do_plot
    title(t, sprintf('Cell %d - raster for %d matched stimulus conditions (%d repeats each)', ...
        cell_id, n_cond, num_reps));
    xlabel(t, 'Time from stimulus onset (s)');
    ylabel(t, 'Trial');
end

end