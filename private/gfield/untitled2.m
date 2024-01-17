%% FIGURE 3

datapath = '/Users/gfield/Analysis/rat/2012-05-02-0/data001/data001';
datapath = '/Users/gfield/Analysis/rat/2012-10-15-0/data000-map/data000-map';
datapath = '/Users/gfield/Analysis/rat/2012-02-23-1/data000/data000';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);

new_purple = [136/256 46/256 114/256];
new_blue = [25/256 101/256 176/256];
new_orange = [230/256 85/256 24/256];


light_purple = [171/256 122/256 150/256];
light_blue = [32/256 164/256 255/256];
light_orange = [237/256 152/256 108/256];

%% Plot ON TCs 

on_type_1_tcs = get_time_courses_matrix(datarun, 'ON type1', 'norm_flag', true);
mean_type_1 = mean(on_type_1_tcs);
on_type_2_tcs = get_time_courses_matrix(datarun, 'ON type2', 'norm_flag', true);
mean_type_2 = mean(on_type_2_tcs);
on_type_3_tcs = get_time_courses_matrix(datarun, 'ON type3', 'norm_flag', true);
mean_type_3 = mean(on_type_3_tcs);

figure(1); clf;
plot([1 30], [0 0], 'k', 'LineWidth', 0.5)
hold on
plot(on_type_1_tcs', 'Color', light_purple, 'LineWidth', 0.25)
plot(mean_type_1, 'Color', new_purple, 'LineWidth', 3)
plot(on_type_2_tcs', 'Color', light_blue, 'LineWidth', 0.5)
plot(mean_type_2, 'Color', new_blue, 'LineWidth', 3)
plot(on_type_3_tcs', 'Color',light_orange, 'LineWidth', 0.5)
plot(mean_type_3, 'Color', new_orange, 'LineWidth', 3)
hold off
axis tight
axis square

print(1, '~/Desktop/ON-tcs.pdf', '-dpdf')

%% ON ISI
isi_duration = 0.054;
isi_bins = 0.001;
isi_length = isi_duration./isi_bins;
datarun = get_interspikeinterval(datarun, {'ON type1', 'ON type2', 'ON type3'}, 'bin_size', isi_bins, 'duration', isi_duration);
on_type1_inds = get_cell_indices(datarun, 'ON type1');
on_type2_inds = get_cell_indices(datarun, 'ON type2');
on_type3_inds = get_cell_indices(datarun, 'ON type3');

on_type1_ISI = zeros(length(on_type1_inds), isi_length+1);
for rgc = 1:length(on_type1_inds)
    tmp_isi = datarun.interspikeinterval{on_type1_inds(rgc)}.probabilities;
    on_type1_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end
on_type2_ISI = zeros(length(on_type2_inds), isi_length+1);
for rgc = 1:length(on_type2_inds)
    tmp_isi = datarun.interspikeinterval{on_type2_inds(rgc)}.probabilities;
    on_type2_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end
on_type3_ISI = zeros(length(on_type3_inds), isi_length+1);
for rgc = 1:length(on_type3_inds)
    tmp_isi = datarun.interspikeinterval{on_type3_inds(rgc)}.probabilities;
    on_type3_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end

tmp_bins = datarun.interspikeinterval{on_type1_inds(1)}.bins;

figure(1); clf;
plot(tmp_bins, on_type1_ISI', 'Color', light_purple)
hold on
plot(tmp_bins, on_type2_ISI', 'Color', light_blue)
plot(tmp_bins, on_type3_ISI', 'Color', light_orange)
plot(tmp_bins, mean(on_type1_ISI)./norm(mean(on_type1_ISI)), 'Color', new_purple, 'LineWidth', 3)
plot(tmp_bins, mean(on_type2_ISI)./norm(mean(on_type2_ISI)), 'Color', new_blue, 'LineWidth', 3)
plot(tmp_bins, mean(on_type3_ISI)./norm(mean(on_type3_ISI)), 'Color', new_orange, 'LineWidth', 3)
axis ([0 0.05 0 0.9])
axis square
axis tight
hold off

print(1, '~/Desktop/ON-ISI.pdf', '-dpdf')

%% OFF ISI
isi_duration = 0.054;
isi_bins = 0.001;
isi_length = isi_duration./isi_bins;
datarun = get_interspikeinterval(datarun, {'OFF type1', 'OFF type2', 'OFF type3'}, 'bin_size', isi_bins, 'duration', isi_duration);
off_type1_inds = get_cell_indices(datarun, 'OFF type1');
off_type2_inds = get_cell_indices(datarun, 'OFF type2');
off_type3_inds = get_cell_indices(datarun, 'OFF type3');

off_type1_ISI = zeros(length(off_type1_inds), isi_length+1);
for rgc = 1:length(off_type1_inds)
    tmp_isi = datarun.interspikeinterval{off_type1_inds(rgc)}.probabilities;
    off_type1_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end
off_type2_ISI = zeros(length(off_type2_inds), isi_length+1);
for rgc = 1:length(off_type2_inds)
    tmp_isi = datarun.interspikeinterval{off_type2_inds(rgc)}.probabilities;
    off_type2_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end
off_type3_ISI = zeros(length(off_type3_inds), isi_length+1);
for rgc = 1:length(off_type3_inds)
    tmp_isi = datarun.interspikeinterval{off_type3_inds(rgc)}.probabilities;
    off_type3_ISI(rgc,:) = tmp_isi./norm(tmp_isi);
end

tmp_bins = datarun.interspikeinterval{on_type1_inds(1)}.bins;

figure(1); clf;
plot(tmp_bins, off_type1_ISI', 'Color', light_purple)
hold on
plot(tmp_bins, off_type2_ISI', 'Color', light_blue)
plot(tmp_bins, off_type3_ISI', 'Color', light_orange)
plot(tmp_bins, mean(off_type1_ISI)./norm(mean(off_type1_ISI)), 'Color', new_purple, 'LineWidth', 3)
plot(tmp_bins, mean(off_type2_ISI)./norm(mean(off_type2_ISI)), 'Color', new_blue, 'LineWidth', 3)
plot(tmp_bins, mean(off_type3_ISI)./norm(mean(off_type3_ISI)), 'Color', new_orange, 'LineWidth', 3)
axis ([0 0.05 0 0.9])
axis square
axis tight
hold off

print(1, '~/Desktop/OFF-ISI.pdf', '-dpdf')

%% Plot OFF TCs 

off_type_1_tcs = get_time_courses_matrix(datarun, 'OFF type1', 'norm_flag', true);
mean_type_1 = mean(off_type_1_tcs);
off_type_2_tcs = get_time_courses_matrix(datarun, 'OFF type2', 'norm_flag', true);
mean_type_2 = mean(off_type_2_tcs);
off_type_3_tcs = get_time_courses_matrix(datarun, 'OFF type3', 'norm_flag', true);
mean_type_3 = mean(off_type_3_tcs);

figure(2); clf;
plot([1 30], [0 0], 'k', 'LineWidth', 0.5)
hold on
plot(off_type_1_tcs', 'Color', light_purple, 'LineWidth', 0.25)
plot(mean_type_1, 'Color', new_purple, 'LineWidth', 3)
plot(off_type_2_tcs', 'Color', light_blue, 'LineWidth', 0.5)
plot(mean_type_2, 'Color', new_blue, 'LineWidth', 3)
plot(off_type_3_tcs', 'Color', light_orange, 'LineWidth', 0.5)
plot(mean_type_3, 'Color', new_orange, 'LineWidth', 3)
hold off
axis tight
axis square

print(2, '~/Desktop/OFF-tcs.pdf', '-dpdf')


%% On type 1

% set array corners for 2012-10-15-0
y_length = 20;
x_length = 2*y_length;
start_x = 22;
start_y = 12;
array_corners = [start_x, start_y; start_x+x_length, start_y; start_x+x_length, start_y+y_length; start_x, start_y+y_length; start_x, start_y];

% set array corners for 2012-02-23-1
y_length = 20;
x_length = 2*y_length;
start_x = 20;
start_y = 6;
array_corners = [start_x, start_y; start_x+x_length, start_y; start_x+x_length, start_y+y_length; start_x, start_y+y_length; start_x, start_y];



temp_cell_type = 'ON type1';
contour_height = 0.6065;
filt_rad = 2.0;

figure(1); clf
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type, contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', true);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/on-type1.eps', '-deps', '-painters')
epsclean('~/Desktop/on-type1.eps','groupSoft',true)

%%
temp_cell_type = 'ON type2';
contour_height = 0.6065;
filt_rad = 0.75;

figure(1); clf
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type,contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', false);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/on-type2.eps', '-deps', '-painters')
epsclean('~/Desktop/on-type2.eps','groupSoft',true)



%%

temp_cell_type = 'ON type3';
contour_height = 0.6065;
filt_rad = 0.75;

figure(1); clf;
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type, contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', false);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/on-type3.eps', '-deps', '-painters')
epsclean('~/Desktop/on-type3.eps','groupSoft',true)



%%

temp_cell_type = 'OFF type1';
contour_height = 0.6065;
filt_rad = 1;

figure(1); clf;
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type, contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', false);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/off-type1.eps', '-deps', '-painters')
epsclean('~/Desktop/off-type1.eps','groupSoft',true)



%%

temp_cell_type = 'OFF type2';
contour_height = 0.6065;
filt_rad = 0.75;

figure(1); clf;
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type, contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', false);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/off-type2.eps', '-deps', '-painters')
epsclean('~/Desktop/off-type2.eps','groupSoft',true)


%%

temp_cell_type = 'OFF type3';
contour_height = 0.6065;
filt_rad = 0.75;

figure(1); clf;
datarun = get_rfs_filtered(datarun,temp_cell_type,'verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));

datarun = get_rf_contours(datarun, temp_cell_type, contour_height, 'norm_params', struct('method', 'peak'), 'verbose', 1);
plot_rf_summaries(datarun, temp_cell_type, 'plot_contours', true, 'contour_colors', 'k','contour_fill', 'k', 'foa', 1, 'clear', false);
plot(array_corners(:,1), array_corners(:,2), 'k')

axis([1 80 1 40])
set(gca,'xtick',[],'ytick',[])
print(1,'~/Desktop/off-type3.eps', '-deps', '-painters')
epsclean('~/Desktop/off-type3.eps','groupSoft',true)


