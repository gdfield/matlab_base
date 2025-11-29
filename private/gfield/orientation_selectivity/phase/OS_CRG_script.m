 
% 2017-07-20-0
cd /Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/TBD/OS-cells/KR_data/
load 2017-07-20-0_data007_004.mat
dg_path = '/Users/gfield/Analysis/2017-07-20-0/data007-map/data007-map';
dg_stim_path = '/Users/gfield/Analysis/2017-07-20-0/stimuli/data007.txt';
crg_path = '/Users/gfield/Analysis/2017-07-20-0/data006_KR-map/data006_KR-map';
crg_stim_path = '/Users/gfield/Analysis/2017-07-20-0/stimuli/data006.txt';
WN_path = '/Users/gfield/Analysis/2017-07-20-0/data004_KR/data004_KR';

%2021-12-28-0
cd ~/Desktop
load 2021-12-28-0.mat
dg_path = '/Users/gfield/Analysis/2021-12-28-0/data001/data001';
dg_stim_path = '/Users/gfield/Analysis/2021-12-28-0/stimuli/s01.txt';
crg_path = '/Users/gfield/Analysis/2021-12-28-0/data003/data003';
crg_stim_path = '/Users/gfield/Analysis/2021-12-28-0/stimuli/s03.txt';
WN_path = '/Users/gfield/Analysis/2021-12-28-0/data000/data000';

%2021-10-07-0
cd /Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/TBD/OS-cells/2021-10-07-0
load 2021-10-07-0_YASS_data001.mat
dg_path = '/Users/gfield/Analysis/2021-10-07-0/Yass/data001/data001';
dg_stim_path = '/Users/gfield/Analysis/2021-10-07-0/stimuli/s01.txt';
crg_path = '/Users/gfield/Analysis/2021-10-07-0/Yass/data002/data002';
crg_stim_path = '/Users/gfield/Analysis/2021-10-07-0/stimuli/s02.txt';
WN_path = '/Users/gfield/Analysis/2021-10-07-0/Yass/data000/data000';

% load DG data
datarun_dg = load_data(dg_path);
datarun_dg = load_neurons(datarun_dg);
datarun_dg = load_params(datarun_dg);
datarun_dg = load_ei(datarun_dg, 'all');
datarun_dg.names.stimulus_path = dg_stim_path;
datarun_dg = load_stim(datarun_dg, 'user_defined_trigger_interval',10);

% extract spatial and temporal periods from datarun
spatial_periods = datarun_dg.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun_dg.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun_dg.stimulus.params.DIRECTION);

os_indices = get_cell_indices(datarun_dg, os_cell_list)

% load CRG data
datarun_crg = load_data(crg_path);
datarun_crg = load_neurons(datarun_crg);
datarun_crg = load_params(datarun_crg);
datarun_crg = load_ei(datarun_crg, 'all');
datarun_crg.names.stimulus_path = crg_stim_path;
datarun_crg = load_stim(datarun_crg);


% load WN data
datarun_wn = load_data(WN_path);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_ei(datarun_wn, 'all');
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');


%% look at selected OS cells
num_rgcs = length(os_cell_list);

cell_indices = get_cell_indices(datarun_dg, os_cell_list);

for rgc = 1:num_rgcs
    % get spikes from an RGC and extract tuning curves
    temp_spike_times = datarun_dg.spikes{os_indices(rgc)};
    
    fig_counter = 1;
    for spat_p = 1:num_sps
        for temp_p = 1:num_tps
            temp_sp = spatial_periods(spat_p);
            temp_tp = temp_periods(temp_p);
            fig_title = ['cell id: ', num2str(datarun_dg.cell_ids(cell_indices(rgc))), ' SP: ', num2str(temp_sp), '; TP: ', num2str(temp_tp)];
            [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_dg.stimulus,...
                                                    'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
            plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_dg.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title)  
            fig_counter = fig_counter + 1;
        end
    end

    % try to map to white noise and plot the spatial RF
    [mapped_list, failed_list] = map_ei(datarun_dg, datarun_wn, 'master_cell_type', datarun_dg.cell_ids(rgc), 'corr_threshold', 0.9);
    temp_cell_index = get_cell_indices(datarun_dg, datarun_dg.cell_ids(rgc));
    wn_cell_id = mapped_list{temp_cell_index};

    % handle case where mapping fails
    figure(10); clf;
    if isempty(wn_cell_id) % this will be empty if a match was not found
        warning('rgc failed to map')
    else % if match found, plot the RF
        plot_rf(datarun_wn, wn_cell_id)
    end
   pause
end

%% 2018
os_cell_list

cell_id = 4457;
temp_index = get_cell_indices(datarun_dg, cell_id);

temp_list = map_ei(datarun_dg, datarun_crg, 'master_cell_type', cell_id);
crg_data_cell_id = temp_list{temp_index};
crg_data_index = get_cell_indices(datarun_crg, crg_data_cell_id);

pixel_phases_80 = [0 80/(360/30) 80/(360/60) 80/(360/90) 80/(360/120) 80/(360/150)];
pixel_phases_120 = [0 120/(360/30) 120/(360/60) 120/(360/90) 120/(360/120) 120/(360/150)];
orientations_list = datarun_crg.stimulus.params.ORIENTATION;

num_phases = length(pixel_phases_80);
num_oris = length(datarun_crg.stimulus.params.ORIENTATION);
% choose spatial period
sp = 80;

spike_counter = zeros(num_oris, num_phases);
num_condits = length(datarun_crg.stimulus.combinations);
for ori_counter = 1:num_oris
    figure(ori_counter); clf;
        for phase_counter = 1:num_phases
        
        g_ori = orientations_list(ori_counter);
        g_phase = pixel_phases_80(phase_counter);

        clear COI
        for condit = 1:num_condits

            sp_check = datarun_crg.stimulus.combinations(condit).SPATIAL_PERIOD;
            ori_check = datarun_crg.stimulus.combinations(condit).ORIENTATION;
            phase_check = datarun_crg.stimulus.combinations(condit).SPATIAL_PHASE;

            if sp_check == sp
                if ori_check == g_ori
                    if round(phase_check) == round(g_phase)
                        COI = condit;
                        break
                    end
                end
            end
        end

        temp_trials = find(datarun_crg.stimulus.trial_list == COI);
        temp_trigs = datarun_crg.stimulus.triggers(temp_trials);


        subplot(1, 6, phase_counter)
        spikes_by_trials = get_raster(datarun_crg.spikes{crg_data_index}, temp_trigs, 'stop', 8, 'foa', ori_counter);

        temp_spike_counter = 0;
        for trial = 1:5
            temp_spike_counter = temp_spike_counter + length(spikes_by_trials{trial});
        end
        spike_counter(ori_counter, phase_counter) = temp_spike_counter;       
    end
end

figure
plot(spike_counter')


%% 2017
os_cell_list

cell_id = 4397;
temp_index = get_cell_indices(datarun_dg, cell_id);

%temp_list = map_ei(datarun_dg, datarun_crg, 'master_cell_type', cell_id);
%crg_data_cell_id = temp_list{temp_index};
%crg_data_index = get_cell_indices(datarun_crg, crg_data_cell_id);

spatial_period = 50;

switch spatial_period
    
    case 25
        pixel_phases = [0 25/(360/30) 25/(360/60) 25/(360/90) 25/(360/120) 25/(360/150)];
        
    case 50
        pixel_phases = [0 50/(360/30) 50/(360/60) 50/(360/90) 50/(360/120) 50/(360/150)];
        
    case 100
        pixel_phases = [0 100/(360/30) 100/(360/60) 100/(360/90) 100/(360/120) 100/(360/150)];

    case 200
        pixel_phases = [0 200/(360/30) 200/(360/60) 200/(360/90) 200/(360/120) 200/(360/150)];

    case 400
        pixel_phases = [0 400/(360/30) 400/(360/60) 400/(360/90) 400/(360/120) 400/(360/150)];
end

orientations_list = [0];

num_phases = length(pixel_phases);
num_oris = length(orientations_list);
% choose spatial period


spike_counter = zeros(num_oris, num_phases);
num_condits = length(datarun_crg.stimulus.combinations);
for ori_counter = 1:num_oris
    figure(ori_counter); clf;
        for phase_counter = 1:num_phases
        
        g_ori = orientations_list(ori_counter);
        g_phase = pixel_phases(phase_counter);

        clear COI
        for condit = 1:num_condits

            sp_check = datarun_crg.stimulus.combinations(condit).SPATIAL_PERIOD;
            ori_check = 0;
            phase_check = datarun_crg.stimulus.combinations(condit).SPATIAL_PHASE;

            if sp_check == spatial_period
                if ori_check == g_ori
                    if round(phase_check) == round(g_phase)
                        COI = condit;
                        break
                    end
                end
            end
        end

        temp_trials = find(datarun_crg.stimulus.trial_list == COI);
        temp_trigs = datarun_crg.stimulus.triggers(temp_trials);


        subplot(1, 6, phase_counter)
        length(temp_trigs)
        spikes_by_trials = get_raster(datarun_crg.spikes{temp_index}, temp_trigs, 'stop', 8, 'foa', ori_counter);

        temp_spike_counter = 0;
        for trial = 1:4
            temp_spike_counter = temp_spike_counter + length(spikes_by_trials{trial});
        end
        spike_counter(ori_counter, phase_counter) = temp_spike_counter;       
    end
end

figure
plot(spike_counter')

