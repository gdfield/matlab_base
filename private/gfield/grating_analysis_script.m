addpath(genpath('/Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code/'))

datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-02-0/data004/data004';

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun.names.stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-02-0/stimuli/s04.txt';

%% process triggers and extract some stim params
trigger_set = round(datarun.triggers);
trig_inds = find(mod(trigger_set, 10) == 0);
user_defined_trigs = datarun.triggers(trig_inds);
datarun = load_stim(datarun, 'user_defined_trigger_set', trig_inds);
num_stim = length(datarun.stimulus.combinations);
grat_dirs = datarun.stimulus.params.DIRECTION;
grat_TPs = datarun.stimulus.params.TEMPORAL_PERIOD;
stim_dur = 8;
num_reps = datarun.stimulus.repetitions;
num_rgcs = length(datarun.cell_ids);



%%
g_sp = datarun.stimulus.params.SPATIAL_PERIOD(1);
for g_dirs = 1:length(datarun.stimulus.params.DIRECTION)
    for g_tp = 1:length(datarun.stimulus.params.TEMPORAL_PERIOD)
        tmp_dir = datarun.stimulus.params.DIRECTION(g_dirs);
        tmp_tp = datarun.stimulus.params.TEMPORAL_PERIOD(g_tp);
        temp_spike_times = get_grating_spike_times(datarun, 'all', stim_dur, 'direction', tmp_dir, 'TP', tmp_tp, 'SP', g_sp);
        gratingrun.direction(g_dirs).temporal_period(g_tp).spike_times = temp_spike_times;
    end
end

%%

[vector_sums_120, vector_mags_120] = get_vector_sums(datarun, 'all', 'TP', 120, 'SP', 240);
[vector_sums_240, vector_mags_240] = get_vector_sums(datarun, 'all', 'TP', 240, 'SP', 240);

scatter((vector_mags_120), (vector_mags_240))




%%
% set x-y cuoff
cutoff_coord = [0.9, 0.25];

x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);


%% Inspect direction tuning

for rgc = 1:length(datarun.cell_ids)
    temp_spikes = datarun.spikes{rgc};
    [tuning_struct, spike_nums, spike_rates] = get_direction_tuning(temp_spikes, datarun.stimulus, 'TP', 60);
    
    datarun.gratings(rgc).tuning_struct = tuning_struct;
    datarun.gratings(rgc).spike_nums = spike_nums;
    datarun.gratings(rgc).spike_rates = spike_rates;

    plot_direction_tuning(tuning_struct, spike_nums, datarun.stimulus,...
            'print', false, 'fig_title', ['cell ', num2str(datarun.cell_ids(rgc))], ...
            'save_path', '~/Desktop/tmp_gratings/', 'save_name', ['cell ', num2str(datarun.cell_ids(rgc))])
end





%% load slave data runs

slave_path_one = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-08-0/data000/data000';
slave_path_two = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-02-0/data001-map/data001-map';
slave_path_wn = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-02-0/data005/data005';

datarun_s1 = load_data(slave_path_one);
datarun_s1 = load_neurons(datarun_s1);
datarun_s1 = load_params(datarun_s1);
datarun_s1 = load_ei(datarun_s1, 'all', 'array_type', 519);
[cell_list_s1, failed_s1] = map_ei(datarun, datarun_s1, 'master_cell_type', ds_cell_ids);

datarun_s2 = load_data(slave_path_two);
datarun_s2 = load_neurons(datarun_s2);
datarun_s2 = load_params(datarun_s2);
datarun_s2 = load_ei(datarun_s2, 'all', 'array_type', 519);
[cell_list_s2, failed_s2] = map_ei(datarun, datarun_s2, 'master_cell_type', ds_cell_ids, 'corr_threshold', 0.9);


datarun_wn = load_data(slave_path_wn);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_ei(datarun_wn, 'all', 'array_type', 519);
[cell_list_wm, failed_wn] = map_ei(datarun, datarun_wn, 'master_cell_type', ds_cell_ids);

num_mapped = length(ds_cell_ids) - length(failed_s2);
master_cell_ids = zeros(num_mapped,1);
slave_id = zeros(num_mapped,1);
num_rgcs = length(cell_list_s2);
cntr = 1;
for rgc = 1:num_rgcs
    if ~isempty(cell_list_s2{rgc})
        master_cell_ids(cntr) = datarun.cell_ids(rgc);
        slave_id(cntr) = cell_list_s2{rgc};
        cntr = cntr +1;
    end
end





