%% Comparison within experiment

datapath = '/Volumes/FieldLab/Chichilnisky/Analysis/2012-02-23-1/kilosort_data000/data000/data000';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);

mark_params.thresh = 3.5;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', mark_params);

%% 
% first cell type
typeA = 'on nc25';
plot_rf_summaries(datarun, typeA, 'plot_fits', true, 'foa', 1)
profiles_A = get_nearest_neighbor_profile(datarun, typeA);

%%
% second cell type
typeB = 'on nc17';
plot_rf_summaries(datarun, typeB, 'plot_fits', true, 'foa', 2)
profiles_B = get_nearest_neighbor_profile(datarun, typeB);

%%
% second cell type
typeC = 'OFF nc81';
plot_rf_summaries(datarun, typeC, 'plot_fits', true, 'foa', 2)
profiles_C = get_nearest_neighbor_profile(datarun, typeC);

%%
figure(5); clf;
plot(profiles_A(1,:)./max(profiles_A(1,:)), 'k')
hold on
plot(profiles_A(2,:)./max(profiles_A(2,:)), 'k')
plot(profiles_B(1,:)./max(profiles_B(1,:)), 'r')
plot(profiles_B(2,:)./max(profiles_B(2,:)), 'r')


plot(profiles_C(1,:)./max(profiles_C(1,:)), 'g')
plot(profiles_C(2,:)./max(profiles_C(2,:)), 'g')

hold off

%% Comparison across experiments

datapath = 'path_to_first_dataset';
datarun1 = load_data(datapath);
datarun1 = load_neurons(datarun1);
datarun1 = load_sta(datarun1, 'load_sta', 'all');
datarun1 = load_params(datarun1);
mark_params.thresh = 3.5;
datarun1 = get_sta_summaries(datarun1, 'all', 'marks_params', mark_params);


datapath = 'path_to_second_dataset';
datarun2 = load_data(datapath);
datarun2 = load_neurons(datarun2);
datarun2 = load_sta(datarun2, 'load_sta', 'all');
datarun2 = load_params(datarun2);
mark_params.thresh = 3.5;
datarun2 = get_sta_summaries(datarun2, 'all', 'marks_params', mark_params);

% first cell type
typeA = 'ON type2';

plot_rf_summaries(datarun1, typeA, 'plot_fits', true, 'foa', 1)
profiles_1 = get_nearest_neighbor_profile(datarun1, typeA);

plot_rf_summaries(datarun2, typeA, 'plot_fits', true, 'foa', 1)
profiles_2 = get_nearest_neighbor_profile(datarun2, typeA);


figure(5); clf;
plot(profiles_1(1,:)./max(profiles_1(1,:)), 'k')
hold on
plot(profiles_1(2,:)./max(profiles_1(2,:)), 'k')
plot(profiles_2(1,:)./max(profiles_2(1,:)), 'r')
plot(profiles_2(2,:)./max(profiles_2(2,:)), 'r')









