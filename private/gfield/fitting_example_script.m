datapath = '/Users/gfield/Analysis/2017-01-16-0/data006_KR/data006_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');


mark_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', mark_params);


cell_indices = get_cell_indices(datarun, 'on brisk transient');

cell_id = datarun.cell_ids(cell_indices(1)); % first off BT RGC

tmp_sta = datarun.stas.stas{cell_indices(1)};


fit_info = fit_sta_sequence(tmp_sta, 'verbose', true);