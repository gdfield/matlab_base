datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-05-19-0/YASS/data002/data002';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

cell_type = 'ON brisk transient';
temp_indices = get_cell_indices(datarun, cell_type);

temp_sta = datarun.stas.stas{temp_indices(1)};

temp_fit_info = fit_sta_sequence(temp_sta, 'verbose', false)



