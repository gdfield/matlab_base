function [tuning_struct, spike_nums, spike_rates] = get_direction_tuning(spike_times, stim_struct, varargin)
%
% usage: [tuning_struct, spike_nums, spike_rates] = get_direction_tuning(spike_times, stim_struct, varargin)
%
% This function takes a list of spike times and a stimulus structure to
% extract the direction tuning of a neuron. 
%
% arguments:  spike_times   -   vector of spike times
%             stim_struct   -   a structure containing the grating parameters,
%                               usually found in datarun.stimulu
%             varargin      -   struct or list of optional parameters (see below)
%
% outputs:    tuning_struct - a cell array of spike times that is
%                             num_directions x num_repeats
%             spike_nums    - 
%                           
%
% optional parameters, their default values, and what they specify:
%
%
% TP              first TP          gets first TP in stim_struct
% SP              first SP         	gets first SP in stim_struct
% stim_duration     8               duration of stimulus in seconds
% contrast        []                only extracts if specified (may cause error if not specified)
% background      []                only extracts if specified
%
%
% Created: GDF, 2019-08-19
%

% ---- BEGIN FUNCTION ---

% parse varargin
p = inputParser;
p.addParameter('TP', stim_struct.params.TEMPORAL_PERIOD(1), @isnumeric);
p.addParameter('SP', stim_struct.params.SPATIAL_PERIOD(1), @isnumeric);
p.addParameter('stim_duration', 8, @isnumeric);
p.addParameter('contrast', []);
p.addParameter('background', []);
p.parse(varargin{:});

% sort information into convenient variable names.
num_repeats = stim_struct.repetitions;
num_stim = length(stim_struct.combinations);
trig_times = stim_struct.triggers;
num_dirs = length(stim_struct.params.DIRECTION);
stim_duration = p.Results.stim_duration;

% find the numeric index to the stim matching the user's input
%stim_indices = zeros(num_dirs, 1);
tuning_struct = cell(num_dirs, num_repeats);
spike_nums = zeros(length(num_dirs), 1);
spike_rates = spike_nums;
for g_dir = 1:num_dirs
    for gstim = 1:num_stim
        if stim_struct.combinations(gstim).TEMPORAL_PERIOD == p.Results.TP && ...
            stim_struct.combinations(gstim).DIRECTION == stim_struct.params.DIRECTION(g_dir) && ...
            stim_struct.combinations(gstim).SPATIAL_PERIOD == p.Results.SP
        stim_index = gstim;
        break
        end
    end
    temp_trigger_inds = find(stim_struct.trial_list == stim_index);

    if length(temp_trigger_inds) ~= num_repeats
        error('insufficient numbr of trials found to match the specified number of repeats');
    end
    
    % loop over rgcs and stim repeats to extract cell array of spike times
    temp_spike_nums = 0;
    for g_rep = 1:num_repeats
        tmp_spike_times = spike_times(spike_times >= trig_times(temp_trigger_inds(g_rep)) &...
                            spike_times < trig_times(temp_trigger_inds(g_rep)) + stim_duration);
        tmp_spike_times = tmp_spike_times - trig_times(temp_trigger_inds(g_rep));

        tuning_struct{g_dir, g_rep} = tmp_spike_times;
        temp_spike_nums = temp_spike_nums + length(tmp_spike_times);
    end
    spike_nums(g_dir) = temp_spike_nums;
    spike_rates(g_dir) = spike_nums(g_dir) ./ num_repeats ./ stim_duration;
end



