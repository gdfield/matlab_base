%% The script is designed to go through datasets listed in the os_datasets 
%% function and allow the user to choose cells that they deem as OS

% Choose and OSI threshold
OSI_thresh = 0.30; % threshold of OS cells
DSI_thresh = 0.25; % threshold for DS cells
response_cor_thresh = 0.3; % threshold on mean correlation across trials
ei_corr_threshold = 0.9; % threshold for EI mapping
filt_rad = 0.75;
fast_plot = true;
wn_map = 1; % decided whether (1) or not (0) to map to white noise


data_list = os_datasets('experiment', 'all');

num_datasets = length(data_list)

dset = 10

for dset = 1:num_datasets
    
    % load stuff
    datarun = load_data(data_list(dset).grating_datapath);
    datarun = load_neurons(datarun);
    datarun = load_params(datarun);
    datarun = load_ei(datarun, 'all');

    % check and establish trigger logic
    temp_stimulus_path = data_list(dset).stimulus_path;
    temp_trigger_interval = data_list(dset).trigger_interval;
    datarun = os_trigger_check(datarun, temp_stimulus_path, temp_trigger_interval); % handles abnormal trigger cases
    
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

    [OSIs, DSIs, mean_corr] = find_OS_RGCs(datarun, 'all', 'verbose', false, 'pause', false);

    os_cell_list = [];
    os_maybe_list = [];
    
    for rgc = 1:num_rgcs

        insig_DSI = all(DSIs(rgc,:,:) < DSI_thresh);
        sig_OSIs = any(OSIs(rgc,:,:) > OSI_thresh);
        sig_cors = any(mean_corr(rgc,:,:) > response_cor_thresh);

        if sig_OSIs && insig_DSI && sig_cors
            fig_counter = 1;
            for spat_p = 1:num_sps
                for temp_p = 1:num_tps
                    temp_sp = spatial_periods(spat_p);
                    temp_tp = temp_periods(temp_p);
                    temp_spike_times = datarun.spikes{rgc};
                    fig_title = ['cell id: ', num2str(datarun.cell_ids(rgc)), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
                    [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                            'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                    plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title, 'fast_plot', fast_plot)  
                    fig_counter = fig_counter + 1;
                end
            end
            fprintf('OSIs are %g \n', OSIs(rgc,:,:))
            fprintf('DSIs are %g \n', DSIs(rgc,:,:))
            fprintf('cors are %g \n', mean_corr(rgc,:,:))

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

