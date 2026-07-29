%% basic data run to generate matched cells from WN to a stimulus set,
% load data from WN and map to other stim types (i.e. gratings/natural
% movies) load in specific cell types of interest into first data run
% (master) i.e. ON parasol and then load in second dataset (slave) and
% map_ei. Use mapped indicies to generate rasters for stimuli set based on
% original classifications. 

% Original Data loading change 'all' to cell type 
datarun = load_data('/Volumes/fieldlab/Array-data/sorted/20260728A/kilosort25/data003/data003');
datarun = load_params(datarun);
datarun.cell_types = order_cell_types(load_txt_cell_types('/Users/mcarleton/Desktop/data003.classification.txt'));
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_ei(datarun, 'all');

% To be mapped Data loading 
datarun1 = load_data('/Volumes/fieldlab/Array-data/sorted/20260728A/kilosort25/data002/data002');
%datarun1 = load_params(datarun1);
datarun1 = load_neurons(datarun1);
datarun1 = load_ei(datarun1, 'all');

%% Mapping from WN to Stimulus set
% Removes nans and makes cell array numeric so it can be used as an input
% in map_ei, both can be set to all since you're working with one cell type
[cell_map, failed, out_mat] = map_ei(datarun, datarun1 ,'master_cell_type', 'all','slave_cell_type', 'all');

% 1. Force the output into a standard numeric matrix if it's packed in a cell array
if iscell(cell_map)
    mapped_cell_ids = cell2mat(cell_map);
end

% Generate master matched IDs if needed (i.e. for spatial correlations with
% RF mapping)
master_idx_matched = find(~cellfun(@isempty, cell_map));
master_ids_matched = datarun.cell_ids(master_idx_matched);

%% Chirp stimulus
chirp_path = '/Volumes/fieldlab/Array-data/sorted/20260715A/kilosort25/data007/data007_ChirpStimulus.mat';
datarun1 = load_chirp(datarun1, chirp_path, 'qi_threshold', 1);

responsive_ids = datarun1.cell_ids(datarun1.chirp.responsive);
plot_chirp_psth(datarun1, responsive_ids(1))
plot_chirp_raster(datarun1, responsive_ids(1))

%Look at mapped_cell_ids to chirp stimulus
browse_chirp(datarun1, mapped_cell_ids)

%Look at chirp versus STA
datarun = get_sta_summaries(datarun, master_ids_matched);
browse_chirp_sta(datarun, datarun1, master_ids_matched, mapped_cell_ids, 'tc_width', 0.22, 'chirp_width', 0.62)

%% plot RF of mapped cells
plot_rf_summaries(datarun,master_ids_matched);

% saves center-of-mass x,y coordinates for correlation analysis
% get the corresponding cell indices for these cell IDs
master_indices_matched = get_cell_indices(datarun, master_ids_matched);

% preallocate: one row per cell, columns = [x, y] center coordinates
rf_centers = nan(length(master_indices_matched), 2);

for i = 1:length(master_indices_matched)
    this_fit = datarun.stas.fits{master_indices_matched(i)};

    if isempty(this_fit) || ~isfield(this_fit, 'mean') || isempty(this_fit.mean)
        warning('No RF fit/mean found for cell ID %d - leaving as NaN.', master_ids_matched(i));
        continue
    end

    rf_centers(i, :) = this_fit.mean(:)';   % [x_center, y_center]
end

% plot center points to check results
figure;
scatter(rf_centers(:,1), rf_centers(:,2), 60, 'filled');
xlabel('X center (pixels)');
ylabel('Y center (pixels)');
title('Receptive field centers');
axis equal;   % keeps spatial proportions correct, since these are x/y coordinates

%% 
% [results, datarun1] = get_all_condition_spike_times(datarun1, cell_ids, trigger_interval, varargin)
% results = get_all_condition_spike_times(datarun, cell_ids, trigger_interval, 'stim_dur', 8)

datarun1.names.stimulus_path = '/Volumes/fieldlab/Array-data/sorted/20260715A/kilosort25/data007/data007_ChirpStimulus.mat';
trigger_set = round(datarun1.triggers);
user_defined_trigs = datarun1.triggers(trig_inds);
datarun1 = load_stim(datarun1, 'user_defined_trigger_set', trig_inds);

%[results, datarun1] = get_all_condition_spike_times(datarun1, mapped_cell_ids, 10);



%% Call testraster to ensure trial structure looks correct

% look at since trial data can input different vargain inputs:
%   testraster(datarun, cell_id)
%   testraster(datarun, cell_id, 'direction', 90)
%   testraster(datarun, cell_id, 'direction', 90, 'temporal_period', 256)
[trial_spike_times, firing_rate, bin_centers] = testraster(datarun1, mapped_cell_ids);



%sequence through a list of data to check rasters: 
%   browse_rasters(datarun, cell_ids)
%   browse_rasters(datarun, cell_ids, 'direction', 90)
%   browse_rasters(datarun, cell_ids, 'direction', 90, 'temporal_period', 256)
browse_rasters(datarun1, mapped_cell_ids);

%% Save rasters for analysis in IgorPro8:
%   saverasters(datarun, output_dir, cell_type)
%   saverasters(datarun, output_dir, cell_type, 'cell_selection', [31 46 91])
%   saverasters(datarun, output_dir, cell_type, 'direction', 90)
%   saverasters(datarun, output_dir, cell_type, 'direction', 90, 'temporal_period', 256)
output_dir = '/Users/mcarleton/Desktop/tutorialdat/outputigor';
saverasters(datarun1, output_dir, mapped_cell_ids, 'direction', 90);

%% DS 

[vector_sums_60, vector_mags_60] = get_vector_sums(datarun1, 'all', 'TP', 60, 'SP', 160);
[vector_sums_240, vector_mags_240] = get_vector_sums(datarun1, 'all', 'TP', 240, 'SP', 160);

% set x-y cuoff
cutoff_coord = [1.5, 1.5];

x_finder = find(vector_mags_60 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun1.cell_ids(selected_indices);
ds_cell_indices = get_cell_indices(datarun1, ds_cell_ids);

for rgc = 1:length(ds_cell_ids)
    temp_spikes = datarun1.spikes{ds_cell_indices(rgc)};
    [tuning_struct, spike_nums, spike_rates] = get_direction_tuning(temp_spikes, datarun1.stimulus, 'TP', 60);

    datarun1.gratings(rgc).tuning_struct = tuning_struct;
    datarun1.gratings(rgc).spike_nums = spike_nums;
    datarun1.gratings(rgc).spike_rates = spike_rates;

    plot_direction_tuning(tuning_struct, spike_nums, datarun1.stimulus,...
        'print', false, 'fig_title', ['cell ', num2str(datarun1.cell_ids(rgc))], ...
        'save_path', '~/Desktop/tmp_gratings/', 'save_name', ['cell ', num2str(datarun1.cell_ids(rgc))])
    pause
end
