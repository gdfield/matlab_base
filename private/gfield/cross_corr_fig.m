datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

%% hOS RGC pair

cell1 = 6483;
cell2 = 5431;

% get cell indices and spikes
temp_cell_indices = get_cell_indices(datarun, [cell1 cell2]);
spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};

% compute ccf and plot
temp_options.dt = 0.005;
temp_options.offset = 0.5;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')
axis([-0.5 0.5 8 13])
axis square
print(1, '~/Desktop/ccf.pdf', '-dpdf')



%% vOS RGCs
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params, 'keep_rf_coms', false);


% cell pair of 
type_num = get_cell_type_nums(datarun, 'OFF type1');
temp_cell_indices = get_cell_indices(datarun, {type_num});

temp_cell_indices = get_cell_indices(datarun, [1817 2343]);

spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};
temp_options.dt = 0.005;
temp_options.offset = 0.5;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')
axis square
axis([-0.5 0.5 8 13])
print(1, '~/Desktop/ccf.pdf', '-dpdf')




%% ON brisk sustained RGCs

cell1 = 4;
cell2 = 6542;

% get cell indices and spikes
temp_cell_indices = get_cell_indices(datarun, [cell1 cell2]);
spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};

% compute ccf and plot
temp_options.dt = 0.005;
temp_options.offset = 0.5;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')

axis([-0.5 0.5 21 38])
axis square
print(1, '~/Desktop/ccf.pdf', '-dpdf')

% get sta summaries and filter
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, [cell1, cell2], 'keep_rf_coms', false, 'marks_params', marks_params);
filt_params.radius = 1;
datarun = get_rfs_filtered(datarun, [cell1, cell2], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% make rf plot for cell1
temp_rf = datarun.stas.filt_rfs{temp_cell_indices(1)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
drawnow 
print(1, '~/Desktop/rf1.pdf', '-dpdf')

% make rf plot for cell2
temp_rf = datarun.stas.filt_rfs{temp_cell_indices(2)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
drawnow 
print(1, '~/Desktop/rf1.pdf', '-dpdf')

%% ON and OFF pair

cell1 = 1398;
cell2 = 1581;

% get cell indices and spikes
temp_cell_indices = get_cell_indices(datarun, [cell1 cell2]);
spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};

% compute ccf and plot
temp_options.dt = 0.005;
temp_options.offset = 0.5;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')
%axis([-0.5 0.5 8 13])
axis square
print(1, '~/Desktop/ccf.pdf', '-dpdf')

% get sta summaries and filter
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, [cell1, cell2], 'keep_rf_coms', false, 'marks_params', marks_params);
filt_params.radius = 1;
datarun = get_rfs_filtered(datarun, [cell1, cell2], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% make rf plot for cell1
temp_rf = datarun.stas.filt_rfs{temp_cell_indices(1)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
drawnow 
print(1, '~/Desktop/rf1.pdf', '-dpdf')

% make rf plot for cell2
temp_rf = -1*datarun.stas.filt_rfs{temp_cell_indices(2)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0 1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
drawnow 
print(1, '~/Desktop/rf1.pdf', '-dpdf')

%% distance ON type 1 cells

cell1 = 4;
cell2 = 5073;

% get cell indices and spikes
temp_cell_indices = get_cell_indices(datarun, [cell1 cell2]);
spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};

% compute ccf and plot
temp_options.dt = 0.001;
temp_options.offset = 0.1;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')
    
axis([-0.5 0.5 21 38])
axis square
print(1, '~/Desktop/ccf.pdf', '-dpdf')

% get sta summaries and filter
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, [cell1, cell2], 'keep_rf_coms', false, 'marks_params', marks_params);
filt_params.radius = 1;
datarun = get_rfs_filtered(datarun, [cell1, cell2], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');


% make rf plot for cell2
temp_rf = datarun.stas.filt_rfs{temp_cell_indices(2)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0 1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
drawnow 
print(1, '~/Desktop/rf1.pdf', '-dpdf')

%%


cell1 = 889;
cell2 = 1191;

% get cell indices and spikes
temp_cell_indices = get_cell_indices(datarun, [cell1 cell2]);
spikes1 = datarun.spikes{temp_cell_indices(1)};
spikes2 = datarun.spikes{temp_cell_indices(2)};

% compute ccf and plot
temp_options.dt = 0.003;
temp_options.offset = 0.3;
[ccf, ccf_time] = compute_ccf_fix(spikes1, spikes2, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')


