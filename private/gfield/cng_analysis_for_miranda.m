%% Control data set (3M)
datapath_wt = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-08-0/data006/data006';
datarun_wt = load_data(datapath_wt);
datarun_wt = load_neurons(datarun_wt);
datarun_wt = load_params(datarun_wt);
datarun_wt = load_sta(datarun_wt, 'load_sta', 'all');

marks_params.thresh = 4.0;
datarun_wt = get_sta_summaries(datarun_wt, 'all', 'marks_params', marks_params);
figure(1); clf;
plot_rf_summaries(datarun_wt, 'ON brisk transient', 'plot_fits', true)
print(1, '~/Desktop/wt_mosaic.pdf', '-dpdf')

wt_tcs = get_time_courses_matrix(datarun_wt, 'ON brisk transient', 'norm_flag', true);
figure(2); clf;
plot(wt_tcs', 'k')
axis square
axis tight
print(2, '~/Desktop/wt_tcs.pdf', '-dpdf')

%% degeneration data set (3M)
datapath_degen = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-04-18-0/data002/data002';
datarun_deg = load_data(datapath_degen);
datarun_deg = load_neurons(datarun_deg);
datarun_deg = load_params(datarun_deg);
datarun_deg = load_sta(datarun_deg, 'load_sta', 'all');

datarun_deg = get_sta_summaries(datarun_deg, 'all', 'marks_params', marks_params);
figure(1); clf;
plot_rf_summaries(datarun_deg, 'ON brisk transient', 'plot_fits', true)
print(1, '~/Desktop/degen_mosaic.pdf', '-dpdf')

degen_tcs = get_time_courses_matrix(datarun_deg, 'ON brisk transient', 'norm_flag', true);
figure(2); clf;
plot(degen_tcs', 'b')
axis square
axis tight
print(2, '~/Desktop/degen_tcs.pdf', '-dpdf')




%%

