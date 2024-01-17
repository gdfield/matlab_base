datapath = '/Users/gfield/Analysis/Flow/2019-01-11-0/data000/data000';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


%%
datapath_sc = '/Users/gfield/Analysis/Flow/2019-01-11-0/data003/data003';
datarun_sc = load_data(datapath_sc);
datarun_sc = load_neurons(datarun_sc);
datarun_sc = load_ei(datarun_sc, 'all', 'array_type', 519);

%%
datapath_ph = '/Volumes/lab/Experiments/Array/Analysis/2019-01-11-0/data005/data005';

datarun_ph = load_data(datapath_ph);
datarun_ph = load_neurons(datarun_ph);
datarun_ph = load_ei(datarun_ph, 'all', 'array_type', 519);


[cell_map, failed_cells] = map_ei(datarun_sc, datarun_ph);

save cell_map cell_map