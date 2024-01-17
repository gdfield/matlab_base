%% OS analysis

% Choose and OSI threshold
OSI_thresh = 0.40; % threshold of OS cells
DSI_thresh = 0.2; % threshold for DS cells
corr_threshold = 0.9; % threshold for EI mapping

% define dataset to look at
% 2012-10-31-0 (good example with lots of OS cells)
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-31-0/data000-1800-7200/data002-map/data002-map';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-31-0/stimuli/s02';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-31-0/data000-1800-7200/data000-map/data000-map';

%SC rats
% 2021-07-13-0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/data000/data000';

% 2021-07-13-1
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/data000/data000';

% 2021-07-15-0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/data000/data000';

% 2021-06-09-0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/data002/data002';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/stimuli/s02.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/data000/data000';

% 2021-06-02-0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/data002/data002';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/stimuli/s02.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/data001/data001';

%OS experiments
%2021-09-09-0 NDF 3.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data000/data000';

%2021-09-09-0 NDF 0.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data003/data003';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s03.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data002/data002';

%2021-09-09-0 NDF 0.0, GABA block
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data004/data004';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s04.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data005/data005';

%2021-09-09-0 NDF, GABA block washout
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data006/data006';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s06.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data007/data007';




%Go through notebooks and find experiemnts with rat data.



% load stuff
datarun = load_data(grating_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

% load stimulus information for gratings
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_interval', 10);

% extract spatial and temporal periods from datarun
spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun.stimulus.params.DIRECTION);


%% load white noise data and map to grating responses
% % load stuff
datarun_wn = load_data(wn_datapath);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
datarun_wn = load_ei(datarun_wn, 'all');

%[mapped_list, failed_list] = map_ei(datarun, datarun_wn);

%% Loop through RGCs and find cells with large OSIs, but small DSIs

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
            [mapped_list, failed_list] = map_ei(datarun, datarun_wn, 'master_cell_type', datarun.cell_ids(rgc), 'corr_threshold', corr_threshold);
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


% save the list of OS cells
cd '~/Desktop/'
save '2021-09-09-0 data003' os_cell_list 


