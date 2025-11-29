function [datarun_dg_h, datarun_dg_low, datarun_wn_high, datarun_wn_low] = load_os_adaptation_data(data_list);

%% load data from vision files

% GRATINGS
% high
datarun_dg_h = load_data(data_list.high_grating_datapath);
datarun_dg_h = load_neurons(datarun_dg_h);
datarun_dg_h = load_params(datarun_dg_h);
datarun_dg_h = load_ei(datarun_dg_h, 'all');

% low
datarun_dg_l = load_data(data_list.low_grating_datapath);
datarun_dg_l = load_neurons(datarun_dg_l);
datarun_dg_l = load_params(datarun_dg_l);
datarun_dg_l = load_ei(datarun_dg_l, 'all');

% NOISE
% high
datarun_wn_h = load_data(data_list.high_wn_datapath);
datarun_wn_h = load_neurons(datarun_wn_h);
datarun_wn_h = load_params(datarun_wn_h);
datarun_wn_h = load_sta(datarun_wn_h, 'load_sta', 'all');
datarun_wn_h = load_ei(datarun_wn_h, 'all');

% low
datarun_wn_l = load_data(data_list.low_wn_datapath);
datarun_wn_l = load_neurons(datarun_wn_l);
datarun_wn_l = load_params(datarun_wn_l);
datarun_wn_l = load_sta(datarun_wn_l, 'load_sta', 'all');
datarun_wn_l = load_ei(datarun_wn_l, 'all');

%% get the grating trigger intervals and parse for each set

% load stimulus information for gratings
datarun_dg_h.names.stimulus_path = data_list(dset).high_stimulus_file;
datarun_dg_h = load_stim(datarun_dg_h, 'user_defined_trigger_interval', 10);


% load stimulus information for agratings
datarun_dg_l.names.stimulus_path = stimulus_path_low;
datarun_dg_l = load_stim(datarun_dg_l, 'user_defined_trigger_set', all_trig_indices_sorted);

dg_mapped_list = map_ei(datarun_dg_h, datarun_dg_l, 'master_cell_type', os_cell_list);


%% map to high light level gratings to high light level WN




% If for 2021-09-09-0 low light level
set_trig_one = find(diff(datarun_dg_l.triggers) > 1.8 & diff(datarun_dg_l.triggers) < 2.3);
set_trig_two = find(diff(datarun_dg_l.triggers) > 4);
set_trig_three = find(diff(datarun_dg_l.triggers) > 1.01 & diff(datarun_dg_l.triggers) < 1.1);
all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
all_trig_indices_sorted = sort(all_trig_indices, 'ascend');

%For 2021-09-23-0 low light level
% all_trig_indices_sorted = [1, find(diff(datarun_dg_l.triggers) > 2.01)'];




% extract spatial and temporal periods from datarun
spatial_periods = datarun_dg_h.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun_dg_h.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun_dg_h.stimulus.params.DIRECTION);


