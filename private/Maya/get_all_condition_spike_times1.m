function [results, datarun1] = get_all_condition_spike_times1(datarun1, cell_ids, trigger_interval, varargin)
%GET_ALL_CONDITION_SPIKE_TIMES  Find stimulus triggers and collect
%spike times for every stimulus combination, for any stimulus class -
%drifting gratings, light steps, contrast response, or anything else.
%
%   results = get_all_condition_spike_times(datarun, cell_ids, trigger_interval)
%   results = get_all_condition_spike_times(datarun, cell_ids, trigger_interval, 'stim_dur', 8)
%   [results, datarun] = get_all_condition_spike_times(datarun, cell_ids, trigger_interval, 'stimulus_path', s_path)
%
% Assumes datarun already has load_data/load_neurons/load_params done
% (i.e. datarun.triggers and datarun.spikes already exist). This
% function finds the triggers marking the start of each stimulus
% presentation, loads the stimulus combinations, and directly collects
% each cell's spike times for each trial of each combination - no
% separate per-stimulus-type spike-extraction function required, since
% the windowing logic (spikes within [trigger, trigger+stim_dur)) is the
% same regardless of what the stimulus parameters are.
%
% REQUIRED INPUTS
%   datarun            - datarun struct with triggers/neurons/params
%                         already loaded
%   cell_ids            - vector of cell IDs to collect spike times for
%   trigger_interval    - time (s) between the start of successive
%                         stimulus presentations (stim duration + any
%                         inter-stimulus gray). For drifting gratings in
%                         the tutorial this was 8 + 2 = 10.
%
% OPTIONAL NAME-VALUE INPUTS
%   'stim_dur'       - duration (s) of the stimulus itself, used to
%                       window spike collection. Defaults to
%                       trigger_interval (i.e. assumes no gray period).
%   'stimulus_path'  - path to the stimulus file. If datarun.stimulus is
%                       already loaded (e.g. you called load_stim
%                       yourself already), omit this and stimulus
%                       loading is skipped.
%
% OUTPUT
%   results  - struct array, one entry per stimulus combination:
%                .params      - struct of that combination's parameter
%                                values (whatever fields your stimulus
%                                has - DIRECTION/TEMPORAL_PERIOD,
%                                CONTRAST, or none at all)
%                .spike_times - cell array [n_cells x n_trials], each
%                                entry a vector of spike times (s)
%                                relative to that trial's onset
%   datarun  - datarun with datarun.stimulus populated (if it wasn't
%              already) and datarun.stimulus.user_stim_dur set

%% parse inputs
p = inputParser;
addParameter(p, 'stim_dur', trigger_interval, @isnumeric);
addParameter(p, 'stimulus_path', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
stim_dur = p.Results.stim_dur;
stimulus_path = p.Results.stimulus_path;

%% find triggers marking the start of each stimulus presentation, and
%% load the stimulus combinations, if not already done

    trigger_set = round(datarun1.triggers);
    trig_inds = find(mod(trigger_set, trigger_interval) == 0);


    datarun1 = load_stim(datarun1, 'user_defined_trigger_set', trig_inds);


datarun1.stimulus.user_stim_dur = stim_dur;



%% collect spikes for every cell, every trial, every stimulus combination
combos = datarun1.stimulus.combinations;
num_stim = length(combos);
num_cells = length(cell_ids);

results = struct('params', {}, 'spike_times', {});

for k = 1:num_stim
    this_combo = combos(k);

    trial_inds = find(datarun1.stimulus.trial_list == k);
    trial_start_times = datarun1.stimulus.triggers(trial_inds);
    n_trials = length(trial_inds);

    spike_times_block = cell(num_cells, n_trials);

    for c = 1:num_cells
        cell_idx = get_cell_indices(datarun1, cell_ids(c));
        cell_spike_times = datarun1.spikes{cell_idx};

        for rep = 1:n_trials
            t0 = trial_start_times(rep);
            t1 = t0 + stim_dur;
            spike_times_block{c, rep} = cell_spike_times(cell_spike_times >= t0 & cell_spike_times < t1) - t0;
        end
    end

    results(k).params = this_combo;
    results(k).spike_times = spike_times_block;
end

fprintf('Collected spike times for %d cell(s) x %d stimulus combination(s).\n', num_cells, num_stim);
end


