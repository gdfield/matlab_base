

datarun{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-06-0/data001/data001');
datarun{1} = load_neurons(datarun{1});
datarun{1} = load_params(datarun{1});
datarun{1} = load_sta(datarun{1}, 'load_sta', 'all');

datarun{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-06-0/data009/data009');
datarun{2} = load_neurons(datarun{2});
datarun{2} = load_params(datarun{2});
datarun{2} = load_sta(datarun{2}, 'load_sta', 'all');

%marks_params.thresh = 4.5;
%datarun{1} = get_sta_summaries(datarun{1}, 'all', 'marks_params', marks_params, 'keep_rf_coms', false);

%% pair 1
% low light level pair
cell1 = 797;
cell2 = 1066;

%high light level pair
cell3 = 796;
cell4 = 1066;


%% pair 2
% low light level pair
cell1 = 2506;
cell2 = 3348;

%high light level pair
cell3 = 2506;
cell4 = 3346;

%%
% %% pair 2
% % low light level pair
% cell1 = 4789;
% cell2 = 5506;
% 
% %high light level pair
% cell3 = 4788;
% cell4 = 5507;
% 
% %% pair 3
% % low light level pair
% cell1 = 5641;
% cell2 = 6766;
% 
% %high light level pair
% cell3 = 5641;
% cell4 = 6586;

%%

cell_indices_low = get_cell_indices(datarun{1}, [cell1 cell2]);
cell_indices_high = get_cell_indices(datarun{2}, [cell3 cell4]);


spikes_one = datarun{1}.spikes{cell_indices_low(1)};
spikes_two = datarun{1}.spikes{cell_indices_low(2)};

spikes_three = datarun{2}.spikes{cell_indices_high(1)};
spikes_four = datarun{2}.spikes{cell_indices_high(2)};

% compute ccf and plot
temp_options.dt = 0.005;
temp_options.offset = 0.25;
[ccf_low, ccf_time_low] = compute_ccf_fix(spikes_one, spikes_two, temp_options);
[ccf_high, ccf_time_high] = compute_ccf_fix(spikes_three, spikes_four, temp_options);

figure(1)
clf
plot(ccf_time_low, ccf_low,'k', 'LineWidth', 3)
hold on
plot(ccf_time_high, ccf_high,'r', 'LineWidth', 3)
hold off
axis square
axis([-0.25 0.25 0 30])
print(1, '~/Desktop/ccf.pdf', '-dpdf')

%%

% get sta summaries and filter
marks_params.thresh = 4.0;
datarun{1} = get_sta_summaries(datarun{1}, [cell1, cell2], 'keep_rf_coms', false, 'marks_params', marks_params);
datarun{2} = get_sta_summaries(datarun{2}, [cell3, cell4], 'keep_rf_coms', false, 'marks_params', marks_params);

filt_params.radius = 0.75;
datarun{1} = get_rfs_filtered(datarun{1}, [cell1, cell2], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
datarun{2} = get_rfs_filtered(datarun{2}, [cell3, cell4], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
indices_low = get_cell_indices(datarun{1}, [cell1, cell2]);
indices_high = get_cell_indices(datarun{2}, [cell4, cell3]);


polarity = -1;

% make rf plot for cell1
figure(1); clf
subplot(2,2,1)
temp_rf = polarity*datarun{1}.stas.filt_rfs{indices_low(1)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
plot_rf_summaries(datarun{1}, [cell1 cell2], 'plot_fit', true, 'scale', 8, 'clear', false)
drawnow 

% make rf plot for cell2
subplot(2,2,2)
temp_rf = polarity*datarun{1}.stas.filt_rfs{indices_low(2)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
plot_rf_summaries(datarun{1}, [cell1 cell2], 'plot_fit', true, 'scale', 8, 'clear', false)
drawnow 


% make rf plot for cell2
subplot(2,2,3)
temp_rf = polarity*datarun{2}.stas.filt_rfs{indices_high(1)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
plot_rf_summaries(datarun{2}, [cell3 cell4], 'plot_fit', true, 'scale', 8, 'clear', false)
drawnow 

% make rf plot for cell2
subplot(2,2,4)
temp_rf = polarity*datarun{2}.stas.filt_rfs{indices_high(2)};
norm_rf = norm_image(temp_rf);
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
axis tight
axis equal
plot_rf_summaries(datarun{2}, [cell3 cell4], 'plot_fit', true, 'scale', 8, 'clear', false)
drawnow 

%% Plot mosaics
low_cell_type_num = get_cell_type_nums(datarun{1}, 'ON biphasic');
high_cell_type_num = get_cell_type_nums(datarun{2}, 'ON biphasic');

mosaic_color = [0 0 0];
figure(1); clf;
plot_rf_fit(datarun{1}, {low_cell_type_num}, 'sd_radius', 1.2, 'edge_color', mosaic_color, 'fill', true, 'alpha', 0.20, 'fill_color',mosaic_color)
plot_rf_fit(datarun{1}, [cell1 cell2], 'sd_radius', 1.2, 'edge_color', mosaic_color, 'fill', true, 'alpha', 0.50, 'fill_color',mosaic_color)
plot_rf_fit(datarun{1}, [2506 3348], 'sd_radius', 1.2, 'edge_color', [0 0 1], 'fill', true, 'alpha', 0.50, 'fill_color',[0 0 1])
axis([0 39 0 30])

print(1, '~/Desktop/rod_mosaic.pdf', '-dpdf')

figure(2); clf;
plot_rf_fit(datarun{2}, {high_cell_type_num}, 'edge_color', mosaic_color, 'fill', true, 'alpha', 0.20, 'fill_color',mosaic_color)
plot_rf_fit(datarun{2}, [cell3 cell4], 'edge_color', mosaic_color, 'fill', true, 'alpha', 0.50, 'fill_color',mosaic_color)
plot_rf_fit(datarun{2}, [2506 3346], 'edge_color', [0 0 1], 'fill', true, 'alpha', 0.50, 'fill_color',[0 0 1])
axis([0 39 0 30])

print(2, '~/Desktop/cone_mosaic', '-dpdf')





