%% OS analysis

% Choose and OSI threshold
OSI_thresh = 0.40; % threshold of OS cells
DSI_thresh = 0.25; % threshold for DS cells
response_cor_thresh = 0.3; % threshold on mean correlation across trials
ei_corr_threshold = 0.9; % threshold for EI mapping
filt_rad = 0.75;
fast_plot = true;
wn_map = 1; % decided whether (1) or not (0) to map to white noise


data_list = os_datasets('experiment', 'all');

num_datasets = length(data_list)

dset = num_datasets

for dset = 1:num_datasets
    
    % load stuff
    datarun = load_data(data_list(dset).grating_datapath);
    datarun = load_neurons(datarun);
    datarun = load_params(datarun);
    datarun = load_ei(datarun, 'all');
    
    if strcmp(data_list(dset).grating_datapath,'/Volumes/gdf/2021-09-09-0/data001/data001')
        % handle these fucked up triggers
        set_trig_one = find(diff(datarun.triggers) > 1.8 & diff(datarun.triggers) < 2.3);
        set_trig_two = find(diff(datarun.triggers) > 4);
        set_trig_three = find(diff(datarun.triggers) > 1.01 & diff(datarun.triggers) < 1.1);
        all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
        all_trig_indices_sorted = sort(all_trig_indices, 'ascend');
    
        % load stimulus information for gratings
        datarun.names.stimulus_path = data_list(dset).stimulus_path;
        datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);
    
    elseif strcmp(data_list(dset).grating_datapath,'/Volumes/gdf/2021-09-23-0/data003/data003')
        % handle these fucked up triggers
        all_trig_indices_sorted = [1, find(diff(datarun.triggers) > 2.01)'];
        
        % load stimulus information for gratings
        datarun.names.stimulus_path = data_list(dset).stimulus_path;
        datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);

    elseif strcmp(data_list(dset).grating_datapath, '/Volumes/gdf/2022-12-21-0/Yass/data001/data001')
        for i=1:length(datarun.triggers)-1
            trig_dif(i) = round(datarun.triggers(i+1)-datarun.triggers(i));
        end
        trial_trig = find(trig_dif==2);
        trial_trig = [1,trial_trig+1];
        datarun.names.stimulus_path = data_list(dset).stimulus_path;
        datarun = load_stim_amr(datarun,'user_defined_trigger_set', trial_trig);
    else
        % load stimulus information for gratings
        datarun.names.stimulus_path = data_list(dset).stimulus_path;
        datarun = load_stim(datarun, 'user_defined_trigger_interval', data_list(dset).trigger_interval);
    end
    
    % extract spatial and temporal periods from datarun
    spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
    temp_periods = datarun.stimulus.params.TEMPORAL_PERIOD;
    num_sps = length(spatial_periods);
    num_tps = length(temp_periods);
    num_dirs = length(datarun.stimulus.params.DIRECTION);
    
    
    %% load white noise data and map to grating responses
    % % load stuff
    datarun_wn = load_data(data_list(dset).wn_datapath);
    datarun_wn = load_neurons(datarun_wn);
    datarun_wn = load_params(datarun_wn);
    datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
    datarun_wn = load_ei(datarun_wn, 'all');
    
    
    
    %% Loop through RGCs and find cells with large OSIs, but small DSIs
    
    num_rgcs = length(datarun.cell_ids);
    
    os_cell_list = [];
    os_maybe_list = [];
    for rgc = 1:num_rgcs
        OSI = zeros(num_sps,num_tps);
        DSI = zeros(num_sps, num_tps);
        mean_cor = zeros(num_sps,num_tps);

        % get spikes from an RGC and extract tuning curves
        temp_spike_times = datarun.spikes{rgc};
        
        for spat_p = 1:num_sps
            for temp_p = 1:num_tps
                [tuning_struct, temp_spike_nums, ~] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
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
                
                % compute reliability (in the preferred
                % orientation/direction)
                [temp_dir_num, temp_rep_num] = size(tuning_struct);
                temp_bin_size = 0.1;
                temp_bins = 0:temp_bin_size:8;
                trial_hists = zeros(length(temp_bins)-1, temp_rep_num);
                temp_spk_times = [];
                for rn = 1:temp_rep_num
                    if isempty(tuning_struct{max_ind, rn})
                        trial_hists(:, rn) = 0;
                    else
                        temp_trial_hists = histcounts(tuning_struct{max_ind, rn},temp_bins);
                        trial_hists(:,rn) = temp_trial_hists ./ norm(temp_trial_hists);
                    end
                end
                cor_matrix = trial_hists' * trial_hists;
                diag_mask = ~eye(size(cor_matrix));
                mean_cor(spat_p, temp_p) = mean(cor_matrix(diag_mask));

                % compute OSI and DSI for each spatial and temporal period
                OSI(spat_p, temp_p) = (max_val - orth_val) ./ (max_val + orth_val);
                DSI(spat_p, temp_p) = (ds_max - null_val) ./ (ds_max + null_val);
            end
        end
        
        % determine if any cross significance thresholds
        insig_DSI = all(DSI(:) < DSI_thresh);
        sig_OSIs = any(OSI(:) > OSI_thresh);
        sig_cors = any(mean_cor(:) > response_cor_thresh);
    
        if sig_OSIs && insig_DSI && sig_cors
%            continue
%        else
%            if insig_DSI
                fig_counter = 1;
                for spat_p = 1:num_sps
                    for temp_p = 1:num_tps
                        temp_sp = spatial_periods(spat_p);
                        temp_tp = temp_periods(temp_p);
                        fig_title = ['cell id: ', num2str(datarun.cell_ids(rgc)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                        [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                        plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title, 'fast', fast_plot)  
                        fig_counter = fig_counter + 1;
                    end
                end
                fprintf('OSIs are %g \n', OSI)
                fprintf('DSIs are %g \n', DSI)
                fprintf('cors are %g \n', mean_cor)

                if wn_map
                    % try to map to white noise and plot the spatial RF
                    mapped_list = map_ei(datarun, datarun_wn, 'master_cell_type', datarun.cell_ids(rgc), 'corr_threshold', ei_corr_threshold);
                    temp_cell_index = get_cell_indices(datarun, datarun.cell_ids(rgc));
                    wn_cell_id = mapped_list{temp_cell_index};
                    
                    % handle case where mapping fails
                    figure(10); clf;
                    if isempty(wn_cell_id) % this will be empty if a match was not found
                        warning('rgc failed to map')
                    else % if match found, plot the RF
                        %plot_rf(datarun_wn, wn_cell_id)
                        temp_rf = get_rf(datarun_wn, wn_cell_id);
                        temp_params = struct('filt_type','gauss','radius',filt_rad);
                        [filt_rf, ~] = rf_filtered(temp_rf, temp_params);
                        norm_rf = norm_image(filt_rf);
                        imagesc(squeeze(norm_rf(:,:,1)))
                        colormap(brewermap([],'RdBu'))
                        caxis([0,1])
                        axis equal  
                        title(num2str(wn_cell_id))
                    end
                end
                
                % get user input to keep or reject cells as OS
                input_check = 0;
                while input_check == 0
                    u_reply = input('keep this cell? y/n/m:', 's');
                    if strcmp(u_reply, 'y')
                        os_cell_list = [os_cell_list, datarun.cell_ids(rgc)];
                        input_check = 1;
                    elseif strcmp(u_reply, 'n')
                        input_check = 1;
                    elseif strcmp(u_reply, 'm')
                        os_maybe_list = [os_maybe_list, datarun.cell_ids(rgc)];
                        input_check = 1;
                    else
                        warning('input was not y/n/m \n')
                    end
                end
        end
    
    end

    % save the list of OS cells
    cd '~/Desktop/os_analysis/'
    split_path = regexp(data_list(dset).grating_datapath,'/','split');
    temp_filename = [split_path{4},'_os'];
    save(temp_filename, 'os_cell_list', 'os_maybe_list');

end

