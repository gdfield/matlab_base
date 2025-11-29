% Choose and OSI threshold
OSI_thresh = 0.6; % threshold of OS cells
DSI_thresh = 0.2; % threshold for DS cells
corr_treshold = 0.9; % threshold for EI mapping

% define dataset to look at
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-15-0/Yass/data002/data002';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-15-0/stimuli/s02.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-15-0/Yass/data000/data000';

% load stuff
datarun = load_data(grating_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

% load stimulus informaation for gratings
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun);
%datarun = load_stim_old(datarun);

% extract spatial and temporal periods from datarun
spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun.stimulus.params.DIRECTION);


%% load white noise data and map to grating responses
% load stuff
datarun_wn = load_data(wn_datapath);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
datarun_wn = load_ei(datarun_wn, 'all');

num_rgcs = length(datarun.cell_ids);

os_cell_list = [];

for rgc = 1:num_rgcs
    OSI = zeros(num_sps,num_tps);
    
    % get spikes from an RGC and extract tuning curves
    temp_spike_times = datarun.spikes{rgc};
    
    for spat_p = 1:num_sps
        for temp_p = 1:num_tps
            [~, temp_spike_nums, ~] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
            
            % compute Orientation-Selective Index
            [ds_max, ds_ind] = max(temp_spike_nums);
            condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs/2, 2]), 2);
            [max_val, max_ind] = max(condensed_tuning);
            
            % handles case of 8 directions
            if num_dirs == 8

                % for OS
                orth_ind = mod(max_ind+2,4);
                if orth_ind == 0
                    orth_ind = 4;
                end
                orth_val = condensed_tuning(orth_ind);
                
                % for DS
                op_ind = mod(ds_ind+4, 8);
                if op_ind == 0
                    op_ind = 1;
                end
                null_val = temp_spike_nums(op_ind);
                
            % handles case of 12 directions    
            elseif num_dirs == 12
                
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
            OSI(spat_p, temp_p) = (max_val - orth_val) ./ (max_val + orth_val);
            DSI(spat_p, temp_p) = (ds_max - null_val) ./ (ds_max + null_val);
        end
    end
    
    % determine if any cross significance thresholds
    insig_DSI = find(DSI(:) < DSI_thresh);
    sig_OSIs = find(OSI(:) > OSI_thresh);

    if isempty(sig_OSIs)
        continue
    else
        if ~isempty(insig_DSI)
            fig_counter = 1;
            for spat_p = 1:num_sps
                for temp_p = 1:num_tps
                    temp_sp = spatial_periods(spat_p);
                    temp_tp = temp_periods(temp_p);
                    fig_title = ['cell id: ', num2str(datarun.cell_ids(rgc)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                    [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                            'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                    plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title)  
                    fig_counter = fig_counter + 1;
                end
            end
            OSI
            
            % try to map to white noise and plot the spatial RF
            [mapped_list, failed_list] = map_ei(datarun, datarun_wn, 'master_cell_type', datarun.cell_ids(rgc), 'corr_threshold', corr_treshold);
            temp_cell_index = get_cell_indices(datarun, datarun.cell_ids(rgc));
            wn_cell_id = mapped_list{temp_cell_index};
            
            % handle case where mapping fails
            figure(10); clf;
            if isempty(wn_cell_id) % this will be empty if a match was not found
                warning('rgc failed to map')
            else % if match found, plot the RF
                plot_rf(datarun_wn, wn_cell_id)
            end
            
            % get user input to keep or reject cells as OS
            input_check = 0;
            while input_check == 0
                u_reply = input('keep this cell? y/n:', 's');
                if strcmp(u_reply, 'y')
                    os_cell_list = [os_cell_list, datarun.cell_ids(rgc)];
                    input_check = 1;
                elseif strcmp(u_reply, 'n')
                    input_check = 1;
                else
                    warning('input was not y/n \n')
                    input_check = 0; % this line isn't needed, but just for bookkeeping
                end
            end
        end
    end
end

