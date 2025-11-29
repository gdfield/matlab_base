%datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-31-0/YASS/data000/data000';
datapath = '/Users/gfield/Analysis/2012-10-31-0/YASS/data000/data000';

datarun = load_data(datapath);

datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.0;

cell_type = {'OFF hOS'};
datarun = get_sta_summaries(datarun, cell_type, 'marks_params', marks_params);

temp_indices = get_cell_indices(datarun, cell_type);

params.radius = 0.75;

cell_list =[];
for rgc = 1:length(temp_indices)
    temp_rf = datarun.stas.rfs{temp_indices(rgc)};

    norm_rf = norm_image(temp_rf);
    
    
    filt_rf = rf_filtered(squeeze(norm_rf(:,:,1)), params);

    imagesc(filt_rf)
    
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    axis equal
    axis tight
    set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
    title(num2str(datarun.cell_ids(temp_indices(rgc))))
    
    user_resp = input('keep this cell y/n', 's');
    if strcmp(user_resp, 'y')
        datarun.cell_ids(temp_indices(rgc))
        cell_list = [cell_list, datarun.cell_ids(temp_indices(rgc))];
    end
end

save 2010-10-31-0-cell_list cell_list


% cd ~/Desktop
% load 2010-10-31-0-cell_list


figure(2)
plot_rf_summaries(datarun, cell_type, 'plot_fits', true)
exportgraphics(gcf, 'hos_mosaic.pdf', 'ContentType', 'vector')

cell_id =472
params.radius = 0.75;
temp_ind = get_cell_indices(datarun, cell_id);
temp_rf = datarun.stas.rfs{temp_ind};
norm_rf = norm_image(temp_rf);
filt_rf = rf_filtered(squeeze(norm_rf(:,:,1)), params);

figure(1)
imagesc(filt_rf)
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
axis equal
axis tight
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
title(num2str(datarun.cell_ids(temp_ind)))
    
%%
datapath = '/Users/gfield/Analysis/2012-10-15-0/YASS/data000/data000';

datarun = load_data(datapath);

datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.0;

cell_type = {'OFF vOS'};
datarun = get_sta_summaries(datarun, cell_type, 'marks_params', marks_params);


