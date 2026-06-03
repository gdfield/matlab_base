% Script to classify different types of OS cells identified as possibly OS
% from the os_cell_finder_script

make_graphics = 0;

% cd to the os analysis directory
cd ~/Desktop/os_analysis/

dset = 1; % this needs to be updated to get that dataset in os_datasets.m

% get the paths to the data
data_list = os_datasets;
path_info = data_list(dset)

% define the data set to inspect
folderName = '2022-12-21-0';
load 2022-12-21-0_os



% make a folder (if it does not already exist) to put plots.
if make_graphics
    if exist(folderName, 'dir') == 7
        disp('Folder is on the path')
    else
        disp(['making folder for ', folderName])
        mkdir(fullfile(pwd, folderName));
    end
end

% cd into this folder.
cd(folderName)

% load stuff grating data
datarun_g = load_data(path_info.grating_datapath);
datarun_g = load_neurons(datarun_g);
datarun_g = load_params(datarun_g);
datarun_g = load_ei(datarun_g, 'all');
% check and establish trigger logic
temp_stimulus_path = data_list(dset).stimulus_path;
temp_trigger_interval = data_list(dset).trigger_interval;
datarun_g = os_trigger_check(datarun_g, temp_stimulus_path, temp_trigger_interval); % handles abnormal trigger cases

%datarun_g.names.stimulus_path = path_info.stimulus_path;
%datarun_g = load_stim(datarun_g, 'user_defined_trigger_interval', path_info.trigger_interval);

% load WN data
datarun_wn = load_data(path_info.wn_datapath);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
datarun_wn = load_ei(datarun_wn, 'all');

% extract spatial and temporal periods from datarun
spatial_periods = datarun_g.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun_g.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun_g.stimulus.params.DIRECTION);
    
% Choose and OSI threshold
%OSI_thresh = 0.40; % threshold of OS cells
%DSI_thresh = 0.25; % threshold for DS cells
response_cor_thresh = 0.3; % threshold on mean correlation across trials
ei_corr_threshold = 0.9; % threshold for EI mapping
filt_rad = 0.75;
fast_plot = true;
wn_map = 1; % decided whether (1) or not (0) to map to white noise


plot_polar_tuning(datarun_g, os_cell_list)

all_possible_os_cells = [os_cell_list, os_maybe_list];

num_rgcs = length(all_possible_os_cells);
os_cell_indices = get_cell_indices(datarun_g, all_possible_os_cells);

[OSIs, DSIs, mean_corr] = find_OS_RGCs(datarun_g, all_possible_os_cells);

for rgc = 1:num_rgcs
   
    temp_spike_times = datarun_g.spikes{os_cell_indices(rgc)};

    fig_counter = 1;
    for spat_p = 1:num_sps
        for temp_p = 1:num_tps
            temp_sp = spatial_periods(spat_p);
            temp_tp = temp_periods(temp_p);
            fig_title = ['cell id-', num2str(os_cell_list(rgc)), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
            [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_g.stimulus,...
                                                    'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
            plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_g.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title, 'fast_plot', true)  
            drawnow
            if make_graphics
                exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
            end
            fig_counter = fig_counter + 1;
        end
    end
    fprintf('OSIs are %g \n', OSIs(rgc,:,:))
    fprintf('DSIs are %g \n', DSIs(rgc,:,:))
    fprintf('cors are %g \n', mean_corr(rgc,:,:))

    if wn_map
        % try to map to white noise and plot the spatial RF
        mapped_list = map_ei(datarun_g, datarun_wn, 'master_cell_type', os_cell_list(rgc), 'corr_threshold', ei_corr_threshold);
        wn_cell_id = mapped_list{os_cell_indices(rgc)};
        
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
            clim([0,1])
            axis equal  
            temp_title = [num2str(os_cell_list(rgc)), '-grating-', num2str(wn_cell_id), '-wn'];
            title(temp_title);
            drawnow
            if make_graphics
                exportgraphics(gcf, [temp_title,'.pdf'], 'ContentType','image')
            end
        end
    end

    hos_simple_cell_list = [];
    vos_complex_cell_list = [];
    hos_thin_cell_list = [];
    vos_thin_cell_list = [];
    unclassified_list = [];
    
% get user input to keep or reject cells as OS
    input_check = 0;
    while input_check == 0
        u_reply = input('identify as type 1, 2, 3, 4, or 0?:');
        if u_reply ==  1
            hos_simple_cell_list = [hos_simple_cell_list, datarun_g.cell_ids(os_cell_indices(rgc))];
            input_check = 1;
        elseif u_reply == 2
            vos_complex_cell_list = [vos_complex_cell_list, datarun_g.cell_ids(os_cell_indices(rgc))];
            input_check = 1;
        elseif u_reply == 3
            hos_thin_cell_list = [hos_thin_cell_list, datarun_g.cell_ids(os_cell_indices(rgc))];
            input_check = 1;
        elseif u_reply == 4
            vos_thin_cell_list = [vos_thin_cell_list, datarun_g.cell_ids(os_cell_indices(rgc))];
            input_check = 1;
        elseif u_reply == 0
            unclassified_list = [unclassified_list, datarun_g.cell_ids(os_cell_indices(rgc))];
            input_check = 1;
    
        else
            warning('input was not supported (0-4)')
        end
    end
end

