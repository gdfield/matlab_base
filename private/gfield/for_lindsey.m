datapath = 'dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-07-18-0/data009-control_KR/data009-control_KR';

datarun = load_data(datapath);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_neurons(datarun);

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


%%

datapath_dg = 'dusom_fieldlab/All_Staff/lab/Experiments/ArrayAnalysis/2018-07-18-0/data009-drug1-map/data009-drug1-map';

datarun_dg = load_data(datapath_dg);
datarun_dg = load_params(datarun_dg);
datarun_dg = load_sta(datarun_dg, 'load_sta', 'all');
datarun_dg = load_neurons(datarun_dg);

marks_params.thresh = 4.0;
datarun_dg = get_sta_summaries(datarun_dg, 'all', 'marks_params', marks_params);


%%

datapath_wash = 'dusom_fieldlab/All_Staff/lab/Experiments/ArrayAnalysis/2018-07-18-0/data009-wash1-map/data009-wash1-map';

datarun_wash = load_data(datapath_wash);
datarun_wash = load_params(datarun_wash);
datarun_wash = load_sta(datarun_wash, 'load_sta', 'all');
datarun_wash = load_neurons(datarun_wash);

marks_params.thresh = 4.0;
datarun_wash = get_sta_summaries(datarun_wash, 'all', 'marks_params', marks_params);
