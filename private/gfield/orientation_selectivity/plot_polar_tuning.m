function plot_polar_tuning(datarun, cell_spec, varargin)

p = inputParser;
p.addParameter('verbose', false, @islogical)
p.parse(varagin{:})



temp_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(temp_indices);

for rgc = 1:num_rgcs

    temp_spike_times = datarun_g.spikes{temp_indices(rgc)};

    [temp_tuning, temp_spike_nums, ~] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));




