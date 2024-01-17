datapath = '/Users/gfield/Analysis/Flow/2018-09-26-0/data007/data007';

datarun = load_data(datapath);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


temp_indices = get_cell_indices(datarun, {'OFF brisk sustained'});
temp_params.radius = 0.5;

for rgc = 1:length(temp_indices)
    temp_rf = -1*get_rf(datarun, datarun.cell_ids(temp_indices(rgc)));
    [temp_rf, params] = rf_filtered(temp_rf, temp_params);

    
    figure(rgc)
    temp_rf = norm_image(temp_rf);
    imagesc(temp_rf(:,:,1))
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    axis equal
    axis tight
    set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
    temp_label = num2str(datarun.cell_ids(temp_indices(rgc)));
    title(temp_label)

end

    
%%

datapath = '/Users/gfield/Analysis/Flow/2018-09-25-0/data000/data000';

datarun = load_data(datapath);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


temp_indices = get_cell_indices(datarun, {'OFF brisk transient'});

for rgc = 1:length(temp_indices)
    temp_rf = -1*get_rf(datarun, datarun.cell_ids(temp_indices(rgc)));

    
    figure(rgc)
    temp_rf = norm_image(temp_rf);
    imagesc(temp_rf(:,:,1))
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    axis equal
    axis tight
    set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
    temp_label = num2str(datarun.cell_ids(temp_indices(rgc)));
    title(temp_label)

end

    