% set server prefix
prefix_flag = 2;
if prefix_flag == 1
    path_prefix = '/Volumes/backup011';
elseif prefix_flag == 2
    path_prefix = '/Users/gfield/Analysis';
elseif prefix_flag == 3
    path_prefix = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis';
end

%2021-09-09-0 NDF 0.0
grating_datapath = '/2021-09-09-0/YASS/data003/data003';
stimulus_path = '/2021-09-09-0/stimuli/s03.txt';
wn_datapath = '/2021-09-09-0/YASS/data002/data002';
load /Users/gfield/Desktop/OS-stuff/OS-cells/2021-09-09-0/2021-09-09-0_data003_yass.mat

%2021-09-23-0 
grating_datapath = '/2021-09-23-0/YASS/data003/data003';
stimulus_path = '/2021-09-23-0/stimuli/s03.txt';
wn_datapath = '/2021-09-23-0/YASS/data004/data004';
load /Users/gfield/Desktop/OS-stuff/OS-cells/2021-09-23-0/2021-09-23-0_data003_yass.mat


%2021-12-28-0 DG and WN
grating_datapath = '/2021-12-28-0/data001/data001';
stimulus_path = '/2021-12-28-0/stimuli/s01.txt';
wn_datapath = '/2021-12-28-0/data000/data000';

%2012-10-31-0
grating_datapath = '/2012-10-31-0/Yass/data002/data002';
stimulus_path = '/2012-10-31-0/stimuli/s02.m';
wn_datapath = '/2012-10-31-0/Yass/data000/data000';
load /Users/gfield/Desktop/OS-stuff/2012-10-31-0-YASS.mat

%2012-10-15-0
grating_datapath = '/2012-10-15-0/Yass/data002/data002';
stimulus_path = '/2012-10-15-0/stimuli/s02.txt';
wn_datapath = '/2012-10-15-0/Yass/data000/data000';
load /Users/gfield/Desktop/OS-stuff/2012-10-15-0-YASS.mat

%2018-11-26-0
grating_datapath = '/2018-11-26-0/data004-map/data004-map';
stimulus_path = '/2018-11-26-0/stimuli/dg.txt';
wn_datapath = '/2018-11-26-0/data003/data003';
load /Users/gfield/Desktop/OS-stuff/OS-cells/KR_data/2018-11-26-0_data003_004.mat


%2019-04-08-0
grating_datapath = '/2019-04-08-0/data002_KR-map/data002_KR-map';
stimulus_path = '/2019-04-08-0/stimuli/dg.txt';
wn_datapath = '/2019-04-08-0/data001_KR/data001_KR';
load /Users/gfield/Desktop/OS-stuff/OS-cells/KR_data/2019-04-08-0_data002.mat


full_grating_datapath = [path_prefix,grating_datapath];
full_stimulus_path = [path_prefix, stimulus_path];
full_wn_path = [path_prefix, wn_datapath];

%%
% load stuff
datarun = load_data(full_grating_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

% load stimulus informaation for gratings
datarun.names.stimulus_path = full_stimulus_path;
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
datarun_wn = load_data(full_wn_path);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
datarun_wn = load_ei(datarun_wn, 'all');
datarun_wn = get_sta_summaries(datarun_wn, 'all');


%% Loop through RGCs and find cells with large OSIs, but small DSIs

num_rgcs = length(os_cell_list);
cell_indices = get_cell_indices(datarun, os_cell_list);

v_os_pow = [];
h_os_pow = [];
for rgc = 1:num_rgcs
    
    % get spikes from an RGC and extract tuning curves
    temp_spike_times = datarun.spikes{cell_indices(rgc)};
    
    fig_counter = 1;
    for spat_p = 1:num_sps
        for temp_p = 1:num_tps
            temp_sp = spatial_periods(spat_p);
            temp_tp = temp_periods(temp_p);
            fig_title = ['cell id: ', num2str(datarun.cell_ids(cell_indices(rgc))), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
            [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                    'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
            plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title)  
            
            tuning_structure{fig_counter} = temp_tuning;           
            fig_counter = fig_counter + 1;
        end
    end   
    
    % try to map to white noise and plot the spatial RF
    [mapped_list, failed_list] = map_ei(datarun, datarun_wn, 'master_cell_type', datarun.cell_ids(rgc), 'corr_threshold', 0.9);
    temp_cell_index = get_cell_indices(datarun, datarun.cell_ids(rgc));
    wn_cell_id = mapped_list{temp_cell_index};                                       

    % handle case where mapping fails
    figure(10); clf;
    if isempty(wn_cell_id) % this will be empty if a match was not found
        warning('rgc failed to map')
    else % if match found, plot the RF
        plot_rf(datarun_wn, wn_cell_id)
    end    

    bin_set = [0:0.001:8];
    spike_counts = zeros(1, length(bin_set)-1);
    %spike_counts = [];
    keeper_flag = input('keep this cell? ', 's');
    if strcmp(keeper_flag, 'y')
        os_flag = input('is this h or v ', 's');
        chosen_orientation = input('choose orientation ');
        chosen_frequency = input('choose plot ');
        spike_set = tuning_structure{chosen_frequency}(chosen_orientation, :);
        num_trials = size(spike_set, 2);
        for trl = 1:num_trials
            temp_counts=histcounts(spike_set{trl}, bin_set);
            spike_counts = spike_counts + temp_counts;
        end
        fft_spikes = fft(spike_counts, [],2);            
        temp_ps = fft_spikes .* conj(fft_spikes);

        figure(11)
        plot(real(temp_ps(1:40)))
        %pause

        start_freq = input('start freq ');
        end_freq = input('end freq ');
        
        freq_pow = sum(temp_ps(start_freq:end_freq)) ./ sum(temp_ps(1:2));

        if strcmp(os_flag, 'v')
            v_os_pow = [v_os_pow, freq_pow]
            save v_os_pow v_os_pow
        elseif strcmp(os_flag, 'h')
            h_os_pow = [h_os_pow, freq_pow]
            save h_os_pow h_os_pow
        else
            error('os input invalid')
        end           
    end  
end




