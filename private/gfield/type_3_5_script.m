% load data, compute RFs and TCs
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);


%2012-10-15-0
offt3 = [62,991,1156,4234,4278,4487,5733,6286,6931];
offt5 = [1246,2253,3695,5116,6260];


TCs = get_time_courses_matrix(datarun, [offt3, offt5]);

[pcomps, pweights] = pca(TCs');

plot(pweights(:,1), pweights(:,2), 'o')

plot(pcomps(:,1), 'k')
hold on
plot(pcomps(:,2), 'r')
hold off

