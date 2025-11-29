% map tuning funcions and spatial RFs across GABA experiments.

local_path = false;

server_prefix = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/';
% local_prefix = '/Users/tierneydaw/Desktop/OS_data_YASS_09-23-2021/';
% 
% if local_path
%     use_prefix = local_prefix;
%     load /Users/tierneydaw/Documents/TDaw/Matlab/OS_cell_lists/2021-09-09-0/2021-09-09-0_data003_yass.mat
% else
    use_prefix = server_prefix;
    load /Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/TBD/OS-cells/2021-09-09-0/2021-09-09-0_data003_yass.mat
% end


%2021-10-07-0 gaba
grating_datapath_control = [use_prefix, 'YASS/data003/data003'];
stimulus_path_control = [use_prefix, 'stimuli/s03.txt'];
grating_datapath_block = [use_prefix, 'YASS/data004/data004'];
stimulus_path_block = [use_prefix, 'stimuli/s04.txt'];

wn_datapath_control = [use_prefix, 'YASS/data002/data002'];
wn_datapath_block = [use_prefix, 'YASS/data005/data005'];


%% map control gratings to GABA block gratings (c = control, b = GABA block)
%load control grating data
datarun_dg_c = load_data(grating_datapath_control);
datarun_dg_c = load_neurons(datarun_dg_c);
datarun_dg_c = load_params(datarun_dg_c);
datarun_dg_c = load_ei(datarun_dg_c, 'all');


% load stimulus information for gratings
datarun_dg_c.names.stimulus_path = stimulus_path_control;
datarun_dg_c = load_stim(datarun_dg_c, 'user_defined_trigger_interval', 10);

% extract spatial and temporal periods from datarun
spatial_periods_c = datarun_dg_c.stimulus.params.SPATIAL_PERIOD;
temp_periods_c = datarun_dg_c.stimulus.params.TEMPORAL_PERIOD;
num_sps_c = length(spatial_periods_c);
num_tps_c = length(temp_periods_c);
num_dirs_c = length(datarun_dg_c.stimulus.params.DIRECTION);

% load GABA block grating data
datarun_dg_b = load_data(grating_datapath_block);
datarun_dg_b = load_neurons(datarun_dg_b);
datarun_dg_b = load_params(datarun_dg_b);
datarun_dg_b = load_ei(datarun_dg_b, 'all');

%For 2021-09-09-0 low light level
% set_trig_one = find(diff(datarun_dg_l.triggers) > 1.8 & diff(datarun_dg_l.triggers) < 2.3);
% set_trig_two = find(diff(datarun_dg_l.triggers) > 4);
% set_trig_three = find(diff(datarun_dg_l.triggers) > 1.01 & diff(datarun_dg_l.triggers) < 1.1);
% all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
% all_trig_indices_sorted = sort(all_trig_indices, 'ascend');

%For 2021-09-23-0 low light level
% all_trig_indices_sorted = [1, find(diff(datarun_dg_l.triggers) > 2.01)'];

% load stimulus information for gratings
datarun_dg_b.names.stimulus_path = stimulus_path_block;
datarun_dg_b = load_stim(datarun_dg_b, 'user_defined_trigger_interval', 10);

% extract spatial and temporal periods from datarun
spatial_periods_b = datarun_dg_b.stimulus.params.SPATIAL_PERIOD;
temp_periods_b = datarun_dg_b.stimulus.params.TEMPORAL_PERIOD;
num_sps_b = length(spatial_periods_b);
num_tps_b = length(temp_periods_b);
num_dirs_b = length(datarun_dg_b.stimulus.params.DIRECTION);


dg_mapped_list = map_ei(datarun_dg_c, datarun_dg_b, 'master_cell_type', os_cell_list);

clear mapped_ids_gratings
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs 
    if ~isempty(dg_mapped_list{rgc}) %EXPLAIN WHATS GOING ON IN THESE 4 LINES
        mapped_ids_gratings(rgc_counter, 1) = datarun_dg_c.cell_ids(rgc);
        mapped_ids_gratings(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_gratings)

%% map control gratings to control WN

datarun_wn_c = load_data(wn_datapath_control);
datarun_wn_c = load_neurons(datarun_wn_c);
datarun_wn_c = load_params(datarun_wn_c);
datarun_wn_c = load_sta(datarun_wn_c, 'load_sta', 'all');
datarun_wn_c = load_ei(datarun_wn_c, 'all');

marks_params.thresh = 4.0;
datarun_wn_c = get_sta_summaries(datarun_wn_c, 'all', 'marks_params', marks_params);

dg_mapped_list = map_ei(datarun_dg_c, datarun_wn_c, 'master_cell_type', os_cell_list);

clear mapped_ids_wn_c
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs
    if ~isempty(dg_mapped_list{rgc})
        mapped_ids_wn_c(rgc_counter, 1) = datarun_dg_c.cell_ids(rgc);
        mapped_ids_wn_c(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_wn_c)


%% map GABA block gratings to GABA block WN

datarun_wn_b = load_data(wn_datapath_block);
datarun_wn_b = load_neurons(datarun_wn_b);
datarun_wn_b = load_params(datarun_wn_b);
datarun_wn_b = load_sta(datarun_wn_b, 'load_sta', 'all');
datarun_wn_b = load_ei(datarun_wn_b, 'all');

marks_params.thresh = 4.0;
datarun_wn_b = get_sta_summaries(datarun_wn_b, 'all', 'marks_params', marks_params);

dg_mapped_list = map_ei(datarun_dg_c, datarun_wn_b, 'master_cell_type', os_cell_list);

clear mapped_ids_wn_b
num_rgcs = length(dg_mapped_list);
rgc_counter = 1;
for rgc = 1:num_rgcs 
    if ~isempty(dg_mapped_list{rgc})
        mapped_ids_wn_b(rgc_counter, 1) = datarun_dg_c.cell_ids(rgc);
        mapped_ids_wn_b(rgc_counter, 2) = dg_mapped_list{rgc};
        rgc_counter = rgc_counter +1;
    end
end

size(mapped_ids_wn_b)


%% identify cells that mapped across all conditions

wn_common = intersect(mapped_ids_wn_b(:,1), mapped_ids_wn_c(:,1));
all_stim_common = intersect(mapped_ids_gratings(:,1), wn_common(:,1));

% get the gratings pair that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_gratings(:,1), all_stim_common); %WHEN IS SOMETHING SANDWHICHED WITH , VS '?
grating_pairs = mapped_ids_gratings(temp_indices,:);;

% get the gratings wn_c that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_wn_c(:,1), all_stim_common);
grating_wn_c_pairs = mapped_ids_wn_c(temp_indices,:);

% get the gratings wn_b that track across all conditions
[~,temp_indices, ~] = intersect(mapped_ids_wn_b(:,1), all_stim_common);
grating_wn_b_pairs = mapped_ids_wn_b(temp_indices,:);



%% plot tuning curves side by side of cells that successfully mapped across drug conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_in_list = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%
control_grating_id =  mapped_ids_gratings(pair_in_list,1)
block_grating_id = mapped_ids_gratings(pair_in_list,2)

control_dg_index = get_cell_indices(datarun_dg_c, control_grating_id);
block_dg_index =  get_cell_indices(datarun_dg_b, block_grating_id);

% get spikes from an RGC and extract tuning curves
temp_spike_times = datarun_dg_c.spikes{control_dg_index}; %HERE, HOW ARE YOU SELECTING FOR RGCS THAT MAP ACROSS ALL CONDITIONS VS RGCS FOUND ONLY IN CONTROL GRATINGS?

% plot tuning curve for cell at control 
fig_counter = 1;
for spat_p = 1:num_sps_c %WHAT ARE THE LEVELS WITHIN THE FOR LOOP? 
    for temp_p = 1:num_tps_c
                temp_sp = spatial_periods_c(spat_p);
                temp_tp = temp_periods_c(temp_p);
                fig_title = ['cell id: ', num2str(datarun_dg_c.cell_ids(control_dg_index)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_c.stimulus,...
                                                        'SP', spatial_periods_c(spat_p), 'TP', temp_periods_c(temp_p));
                plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_c.stimulus,...
                            'fig_num', fig_counter, 'fig_title', fig_title, 'print_for_fig', false)  
                fig_counter = fig_counter + 1;                                                    

        % compute Orientation-Selective Index
        [ds_max, ds_ind] = max(temp_spike_nums); 
        condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs_c/2, 2]), 2); %OK WHAT'S GOING ON HERE?????
        [max_val, max_ind] = max(condensed_tuning);
            
                
            % handles case of 12 directions    
            if num_dirs_c == 12
                
                % for OS
                orth_ind = mod(max_ind+3,6);
                if orth_ind == 0
                    orth_ind = 6;
                end
                orth_val = condensed_tuning(orth_ind);
                
                % for DS
                op_ind = mod(ds_ind+6, 12);
                if op_ind == 0
                    op_ind = 1;
                end
                null_val = temp_spike_nums(op_ind);
                
            end
            
            % compute OSI and DSI for each spatial and temporal period
            OSI_c(spat_p, temp_p) = (max_val - orth_val) ./ (max_val + orth_val);        
   end
end
OSI_c

% get spikes from an RGC and extract tuning curves
temp_spike_times = datarun_dg_b.spikes{block_dg_index};

% extract spatial and temporal periods from datarun
spatial_periods_b = datarun_dg_b.stimulus.params.SPATIAL_PERIOD;
temp_periods_b = datarun_dg_b.stimulus.params.TEMPORAL_PERIOD;
num_sps_b = length(spatial_periods_b);
num_tps_b = length(temp_periods_b);
num_dirs_b = length(datarun_dg_c.stimulus.params.DIRECTION);

% plot tuning curve for cell at low light level
fig_counter = 5;
for spat_p = 1:num_sps_b
    for temp_p = 1:num_tps_b
                temp_sp = spatial_periods_b(spat_p);
                temp_tp = temp_periods_b(temp_p);
                fig_title = ['cell id: ', num2str(datarun_dg_b.cell_ids(block_dg_index)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_b.stimulus,...
                                                        'SP', spatial_periods_b(spat_p), 'TP', temp_periods_b(temp_p));
                plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_b.stimulus,...
                            'fig_num', fig_counter, 'fig_title', fig_title, 'print_for_fig', false)  
                fig_counter = fig_counter + 1;                                                    

        % compute Orientation-Selective Index
        [ds_max, ds_ind] = max(temp_spike_nums);
        condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs_b/2, 2]), 2);
        [max_val, max_ind] = max(condensed_tuning);

               % handles case of 12 directions    
            if num_dirs_b == 12
                
                % for OS
                orth_ind = mod(max_ind+3,6);
                if orth_ind == 0
                    orth_ind = 6;
                end
                orth_val = condensed_tuning(orth_ind);
                
                % for DS
                op_ind = mod(ds_ind+6, 12);
                if op_ind == 0
                    op_ind = 1;
                end
                null_val = temp_spike_nums(op_ind);
                
            end
            
            % compute OSI and DSI for each spatial and temporal period
            OSI_b(spat_p, temp_p) = (max_val - orth_val) ./ (max_val + orth_val);     
    
    end
end
OSI_b

% plot the RF at the control condition
plot_rf(datarun_wn_c, grating_wn_c_pairs(pair_in_list,2), 'foa', 100, 'fit', false)

% temp_rf = datarun_wn_c.stas.rfs{get_cell_indices(datarun_wn_c, grating_wn_c_pairs(pair_in_list,2))};
% expand_rf = repmat(temp_rf, [1 1 3]);
% imagesc(expand_rf)
% colormap(brewermap([], 'RdBu'))

% plot the RF at the GABA block condition
plot_rf(datarun_wn_b, grating_wn_b_pairs(pair_in_list,2), 'foa', 101, 'fit', false)

figure(111)
plot_ei(datarun_dg_c, control_grating_id)
figure(112)
plot_ei(datarun_dg_b, block_grating_id)
figure(113)
plot_ei(datarun_wn_c, grating_wn_c_pairs(pair_in_list,2))
figure(114)
plot_ei(datarun_wn_b, grating_wn_b_pairs(pair_in_list,2))


%% If you like the cell, print shit out to make a figure:

print(100, '~/Desktop/rf_control.pdf', '-dpdf')
print(101, '~/Desktop/rf_block.pdf', '-dpdf')

% choose spatial period and temporal period indices that are best liked
% NOTE on Directions:
% "0" degrees is 3:00, and the angles go counter clockwise

%for hOS cell 6 (2021-09-09-0 NDF0 control v block), used SP = 80 and TP = 60
%for vOS cell 11 (2021-09-09-0 NDF0 control v block), used SP = 80 and TP = 60


SP = 80;
TP = 60; 
temp_sp = find(spatial_periods == SP);
temp_tp = find(temp_periods == TP);
% print for control condition
temp_spike_times = datarun_dg_c.spikes{control_dg_index};
[temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_c.stimulus,...
                                        'SP', spatial_periods(temp_sp), 'TP', temp_periods(temp_tp));
plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_c.stimulus,...
            'fig_num', fig_counter, 'print_for_fig', true)  
    
        
% print for block condition 
temp_spike_times = datarun_dg_b.spikes{block_dg_index};
[temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg_b.stimulus,...
                                        'SP', spatial_periods(temp_sp), 'TP', temp_periods(temp_tp));
plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg_b.stimulus,...
            'fig_num', fig_counter, 'print_for_fig', true)  
        
        




