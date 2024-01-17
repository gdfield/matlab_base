%%

datapath = '/Users/gdfield/Analysis/2013-10-28-0/data001/data001';
datapath = '/Users/gdfield/Analysis/2013-10-30-0/data001/data001';
datapath = '/Users/gdfield/Analysis/2014-10-13-0/data001/data001';
datapath = '/Users/gdfield/Analysis/2014-10-15-0/mapped/data001/unclean/data001/data001';

datapath = '/Users/gdfield/Analysis/2014-10-15-0/mapped/data001/data001';

datapath = '/Users/gdfield/Analysis/2012-10-15-0/data000-map/unclean/data000/data000-map';


%% Off type 5
datapath = '/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map';
marks_params.thresh = 4.5;
temp_params.radius = 0.75;
cell_type = {'OFF unclassified nc5'};

% pseudo_offt5 = [286 1970 2582 3528 4279 5987 7638];
% offt6 = [302 1007 4966];
% cell_type = offt6;
% 
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

datarun = get_sta_summaries(datarun, cell_type, 'marks_params', marks_params);

cell_indices = get_cell_indices(datarun, cell_type);
num_rgcs = length(cell_indices);

% filter the rfs by a Gaussian 
filt_params.radius = 0.75;
datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
datarun = get_rf_contours(datarun, cell_type, 4, 'save_cont_params', 'cont_rfs', 'rfs', 'filt_rfs');

plot_rf_summaries(datarun, cell_type, 'plot_fits', false, 'foa', -1,...
                'bg_color', [], 'clear', true, 'label', true,...
                'plot_contours', true, 'contour_fill', 'r', 'contour_alpha', 0)

%% get the "surround" countours

datarun.stas.filt_rfs2 = cell(length(datarun.cell_ids), 1);
for rgc = 1:num_rgcs
    temp_rf = datarun.stas.filt_rfs{cell_indices(rgc)};
    
    temp_std = robust_std(temp_rf(:));
    temp_thresh = temp_std * 3.5;
    low_thresh = 0 - temp_thresh;
    high_thresh = 0 + temp_thresh;
    threshed_pixel_indices = find(temp_rf > low_thresh & temp_rf < high_thresh);
    temp_rf(threshed_pixel_indices) = 0.0;
      
    
    temp_rf(temp_rf >0) = 0;
    datarun.stas.filt_rfs2{cell_indices(rgc)} = -1*temp_rf;
end
datarun = get_rf_contours(datarun, cell_type, 4, 'save_cont_params', 'cont_rfs', 'rfs', 'filt_rfs2', 'save_name', 'rf_contours2');


plot_rf_summaries(datarun, cell_type, 'plot_fits', false, 'foa', -1,...
                'bg_color', [], 'clear', true, 'label', false,...
                'plot_contours', true, 'contours_field', 'rf_contours2',...
                'contour_colors', 'kb', 'contour_fill', 'b', 'contour_alpha', 0.3)



plot_rf_summaries(datarun, [1731 1816 3288 3678], 'plot_fits', false, 'foa', 2,...
    'bg_color', [], 'clear', false, 'label', true,...
    'plot_contours', true, 'contours_field', 'rf_contours',...
    'contour_alpha', 0.3)

            
%% Plot RF Profiles

summed_rf = zeros(size(datarun.stas.rfs{cell_indices(1)}));
figure(1); clf;
for rgc = 1:num_rgcs
    temp_rf = datarun.stas.rfs{cell_indices(rgc)};
    filt_rf = rf_filtered(temp_rf, temp_params);

    temp_std = robust_std(filt_rf(:));
    temp_thresh = temp_std * 0.5;
    low_thresh = 0 - temp_thresh;
    high_thresh = 0 + temp_thresh;
    threshed_pixel_indices = find(filt_rf > low_thresh & filt_rf < high_thresh);
    filt_rf(threshed_pixel_indices) = 0.0;
       
    summed_rf = filt_rf + summed_rf;
    
    norm_rf = norm_image(filt_rf);
   
    subplot(1,1,1)
    imagesc(matrix_scaled_up(squeeze(norm_rf(:,:,1)),8))
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    axis equal
    axis tight
    xlabel(num2str(datarun.cell_ids(cell_indices(rgc))));
    set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
    pause
end

norm_summed_rf = norm_image(summed_rf);

imagesc(matrix_scaled_up(squeeze(norm_summed_rf(:,:,1)),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
axis equal
axis tight
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

plot_rf_summaries(datarun, cell_type, 'plot_fits', false, 'foa', 2,...
                'bg_color', [], 'clear', true, 'label', false,...
                'plot_contours', true, 'contour_fill', 'r', 'contour_alpha', 0)
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

plot_rf_summaries(datarun, cell_type, 'plot_fits', false, 'foa', 3,...
                'bg_color', [], 'clear', true, 'label', false,...
                'plot_contours', true, 'contours_field', 'rf_contours2',...
                'contour_colors', 'b', 'contour_fill', 'b', 'contour_alpha', 0)
            set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
%%

print(1, '~/Desktop/rfs.pdf', '-dpdf')
print(2, '~/Desktop/center.pdf', '-dpdf')
print(3, '~/Desktop/surround.pdf', '-dpdf')


            
            
