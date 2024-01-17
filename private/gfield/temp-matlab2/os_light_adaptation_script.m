% map tuning funcions and spatial RFs across light levels.


%2021-09-09-0 NDF 0.0
load /Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/TBD/OS-cells/2021-09-09-0/2021-09-09-0-data003-yass.mat

grating_datapath_high = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data003/data003';
stimulus_path_high = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s03.txt';
grating_datapath_low = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data001/data001';
stimulus_path_low = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';

wn_datapath_high = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data002/data002';
wn_datapath_low = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data000/data000';


%% map high light level gratings to low light levels gratings
datarun_dg_h = load_data(grating_datapath_high);
datarun_dg_h = load_neurons(datarun_dg_h);
datarun_dg_h = load_params(datarun_dg_h);
datarun_dg_h = load_ei(datarun_dg_h, 'all');


% load stimulus information for gratings
datarun_dg_h.names.stimulus_path = stimulus_path_high;
datarun_dg_h = load_stim(datarun_dg_h, 'user_defined_trigger_interval', 10);

% extract spatial and temporal periods from datarun
spatial_periods = datarun_dg_h.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun_dg_h.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun_dg_h.stimulus.params.DIRECTION);

% load low light level grating data
datarun_dg_l = load_data(grating_datapath_low);
datarun_dg_l = load_neurons(datarun_dg_l);
datarun_dg_l = load_params(datarun_dg_l);
datarun_dg_l = load_ei(datarun_dg_l, 'all');

set_trig_one = find(diff(datarun_dg_l.triggers) > 1.8 & diff(datarun_dg_l.triggers) < 2.3);
set_trig_two = find(diff(datarun_dg_l.triggers) > 4);
set_trig_three = find(diff(datarun_dg_l.triggers) > 1.01 & diff(datarun_dg_l.triggers) < 1.1);
all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
all_trig_indices_sorted = sort(all_trig_indices, 'ascend');

% load stimulus information for gratings
datarun_dg_l.names.stimulus_path = stimulus_path_low;
datarun_dg_l = load_stim(datarun_dg_l, 'user_defined_trigger_set', all_trig_indices_sorted);

dg_mapped_list = map_ei(datarun_dg_h, datarun_dg_l, 'master_cell_type', os_cell_list);

clear mapped_ids_gratings
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs 
    if ~isempty(dg_mapped_list{rgc})
        mapped_ids_gratings(rgc_counter, 1) = datarun_dg_h.cell_ids(rgc);
        mapped_ids_gratings(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_gratings)

%% map to high light level gratings to high light level WN

datarun_wn_h = load_data(wn_datapath_high);
datarun_wn_h = load_neurons(datarun_wn_h);
datarun_wn_h = load_params(datarun_wn_h);
datarun_wn_h = load_sta(datarun_wn_h, 'load_sta', 'all');
datarun_wn_h = load_ei(datarun_wn_h, 'all');

dg_mapped_list = map_ei(datarun_dg_h, datarun_wn_h, 'master_cell_type', os_cell_list);

clear mapped_ids_wn_h
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs 
    if ~isempty(dg_mapped_list{rgc})
        mapped_ids_wn_h(rgc_counter, 1) = datarun_dg_h.cell_ids(rgc);
        mapped_ids_wn_h(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_wn_h)


%% map to high light level gratings to high light level WN

datarun_wn_l = load_data(wn_datapath_low);
datarun_wn_l = load_neurons(datarun_wn_l);
datarun_wn_l = load_params(datarun_wn_l);
datarun_wn_l = load_sta(datarun_wn_l, 'load_sta', 'all');
datarun_wn_l = load_ei(datarun_wn_l, 'all');

dg_mapped_list = map_ei(datarun_dg_h, datarun_wn_l, 'master_cell_type', os_cell_list);

clear mapped_ids_wn_l
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs 
    if ~isempty(dg_mapped_list{rgc})
        mapped_ids_wn_l(rgc_counter, 1) = datarun_dg_h.cell_ids(rgc);
        mapped_ids_wn_l(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_wn_l)


%% identify cells that mapped across all conditions

wn_common = intersect(mapped_ids_wn_l(:,1), mapped_ids_wn_h(:,1));
all_stim_common = intersect(mapped_ids_gratings(:,1), wn_common(:,1));

% get the gratings pair that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_gratings(:,1), all_stim_common);
grating_pairs = mapped_ids_gratings(temp_indices,:);;

% get the gratings wn_h that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_wn_h(:,1), all_stim_common);
grating_wn_h_pairs = mapped_ids_wn_h(temp_indices,:);

% get the gratings wn_l that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_wn_l(:,1), all_stim_common);
grating_wn_l_pairs = mapped_ids_wn_l(temp_indices,:);



%% plot tuning curves side by side of cells that successfully mapped across
% light levels

%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_in_list = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%
high_grating_id = grating_pairs(pair_in_list,1)
low_grating_id = grating_pairs(pair_in_list,2)

high_dg_index = get_cell_indices(datarun_dg_h, high_grating_id);
low_dg_index =  get_cell_indices(datarun_dg_l, low_grating_id);



% get spikes from an RGC and extract tuning curves
temp_spike_times = datarun_dg_h.spikes{high_dg_index};
   

% plot tuning curve for cell at high light level
fig_counter = 1;
for spat_p = 1:num_sps
    for temp_p = 1:num_tps
                temp_sp = spatial_periods(spat_p);
                temp_tp = temp_periods(temp_p);
                fig_title = ['cell id: ', num2str(datarun_dg_h.cell_ids(high_dg_index)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_h.stimulus,...
                                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_h.stimulus,...
                            'fig_num', fig_counter, 'fig_title', fig_title, 'print_for_fig', false)  
                fig_counter = fig_counter + 1;                                                    

        % compute Orientation-Selective Index
        [ds_max, ds_ind] = max(temp_spike_nums);
        condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs/2, 2]), 2);
        [max_val, max_ind] = max(condensed_tuning);
   end
end

%% get spikes from an RGC and extract tuning curves
temp_spike_times = datarun_dg_l.spikes{low_dg_index};

% plot tuning curve for cell at low light level
fig_counter = 5;
for spat_p = 1:num_sps
    for temp_p = 1:num_tps
                temp_sp = spatial_periods(spat_p);
                temp_tp = temp_periods(temp_p);
                fig_title = ['cell id: ', num2str(datarun_dg_l.cell_ids(low_dg_index)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_l.stimulus,...
                                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_l.stimulus,...
                            'fig_num', fig_counter, 'fig_title', fig_title, 'print_for_fig', false)  
                fig_counter = fig_counter + 1;                                                    

        % compute Orientation-Selective Index
        [ds_max, ds_ind] = max(temp_spike_nums);
        condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs/2, 2]), 2);
        [max_val, max_ind] = max(condensed_tuning);
   end
end


% plot the RF at the high light level
plot_rf(datarun_wn_h, grating_wn_h_pairs(pair_in_list,2), 'foa', 100, 'fit', false)

% plot the RF at the low light level
plot_rf(datarun_wn_l, grating_wn_l_pairs(pair_in_list,2), 'foa', 101, 'fit', false)



%% If you like the cell, print shit out to make a figure:

print(100, '~/Desktop/rf_high.pdf', '-dpdf')
print(101, '~/Desktop/rf_low.pdf', '-dpdf')

% choose spatial period and temporal period indices that are best liked
% NOTE on Directions:
% "0" degrees is 3:00, and the angles go counter clockwise

SP =160;
TP = 60; 
temp_sp = find(spatial_periods == SP);
temp_tp = find(temp_periods == TP);
% print for high light level
temp_spike_times = datarun_dg_h.spikes{high_dg_index};
[temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_h.stimulus,...
                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_h.stimulus,...
            'fig_num', fig_counter, 'print_for_fig', true)  
    
        
% print for low light level
temp_spike_times = datarun_dg_l.spikes{low_dg_index};
[temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_l.stimulus,...
                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_l.stimulus,...
            'fig_num', fig_counter, 'print_for_fig', true)  
        
        




