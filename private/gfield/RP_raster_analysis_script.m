% load data
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2014-06-16-0/data004/data004');
datarun = load_neurons(datarun);
datarun = load_params(datarun);

% choose triggers that identify each white noise repeat.
temp_triggers = datarun.triggers(1:6:1080-5);

% choose which cell
cell_index = 105;

% get and plot rast for cell identified by cell_index
spike_trials = get_raster(datarun.spikes{cell_index}, temp_triggers, 'plot', true, 'tic_color', 'r');





%%


