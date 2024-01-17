%% Info analysis on WN repeats

% define data path and load data
datapath = [server_path, '2019-07-18-0/data004/data004'];
datarun = load_data(datapath);
datarun = load_neurons(datarun);

%% Plot a raster

% extract triggers
diff_trigs = diff(datarun.triggers);
trig_inds = find(diff_trigs > 1.67);

cell_id = 811;
cell_index = get_cell_indices(datarun, cell_id);

temp_spike_times = datarun.spikes{cell_index};
temp_raster = get_raster(temp_spike_times, datarun.triggers(trig_inds));

%% Info calculation

% bin spike train
bin_size = 0.005; % define bin size in seconds
rep_length = 9; % define trial duration in seconds
bg_buffer = 10; % skip this number trials at beginning of recordings

trial_bins = 0:bin_size:rep_length; % edges of bins
num_bins = length(trial_bins)-1; % number of bins (-1 because bins are defined by edges);
tmp_spikes = datarun.spikes{cell_index}; % make spikes train a vector


% bin spikes by time and trial
% initialize binned spikes matrix
binned_spikes = zeros((length(trig_inds) - bg_buffer+1), num_bins); 
for trl = bg_buffer:length(trig_inds)
    tmp_start = datarun.triggers(trig_inds(trl));
    tmp_inds = find(tmp_spikes >= tmp_start &...
                    tmp_spikes < (tmp_start + rep_length));
    epoch_spikes = tmp_spikes(tmp_inds);            
    
    %tracker = tmp_start -            
   tmp_binned_spikes = histcounts(epoch_spikes - tmp_start , [0:bin_size:rep_length]);
   
   binned_spikes(trl-bg_buffer+1,:) = tmp_binned_spikes;
end


word_size = 6;
cell_Info = get_MI_from_resp(binned_spikes, trial_bins, word_size, datarun.triggers(trig_inds(bg_buffer:length(trig_inds)))




