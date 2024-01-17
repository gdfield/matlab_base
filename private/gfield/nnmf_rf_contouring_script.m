% load data, compute RFs and TCs
%2012-10-10-1
datarun = load_data('/Users/gdfield/Analysis/2012-10-10-1/data001-3600-7200s/data001-map/data001-map');
offt3 = [512 1368 2750 4188 4487 5732 6170 6511 6588];
offt5 = [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353,3874];

%2012-10-15-0
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
offt3 = [62,991,1156,4234,4278,4487,5733,6286,6931];
offt5 = [1246,2253,3695,5116,6260];

%2012-10-31
datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map');
offt3 = [1577,1954,2659,7324,6033,4383,2884];
offt5 = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376];


datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

% book keeping
cell_type_num = get_cell_type_nums(datarun, 'on type2');
cell_indices = get_cell_indices(datarun, {cell_type_num});

num_rgcs = length(cell_indices);
for rgc = 1:num_rgcs

    % get the space-time STA
    temp_sta = datarun.stas.stas{cell_indices(rgc)};

    % squeeze and reshape the sta for SVD
    temp_sta = squeeze(temp_sta);
    temp_sta = reshape(temp_sta, [3200,30]);

    % perform SVD
    [U, S, V] = svd(temp_sta);
    
    [W, H] = nnmf(temp_sta,5);
    
    temp_rf = reshape(W(:,1), [40,80]);
    
    imagesc(temp_rf)
    %drawnow
    %pause
    
    datarun.stas.nnmf_rfs{cell_indices(rgc)} = temp_rf;
  
end

filt_params.filt_type = 'gauss';
filt_params.radius = 0.75;

summed_rfs = zeros(40,80);
clear mosaic_contours
for rgc = 1:length(cell_indices)
   
    % filter RF
    temp_rf = datarun.stas.nnmf_rfs{cell_indices(rgc)};
    filt_summary = rf_filtered(temp_rf, filt_params);

    % normalized the filtered summary
    filt_summary = filt_summary ./ max(filt_summary(:));

    % plot filtered RF
     figure(1)
     imagesc(filt_summary); colormap jet
     title(num2str(rgc))
    
    % get contours for RF
    contour_levels = 0.3;
    temp_contour_polygons = rf_contours(double(filt_summary), contour_levels);
    plot_polygon_struct(temp_contour_polygons{1}, 'facecolor', [1 1 1], 'bgcolor', 'k', 'alpha', 0.1, 'linecolor', 'w', 'LineWidth', 1);
    
    pause
    mosaic_contours{rgc} = temp_contour_polygons{1};
    
    % threshold filtered RF
    filt_summary(filt_summary(:) < 0.05) = 0;
    
    if ismember(rgc, [1 2 4 10])
        figure(2)
        summed_rfs = summed_rfs + filt_summary;
        imagesc(summed_rfs); colormap jet
        plot_polygon_struct(temp_contour_polygons{1}, 'facecolor', [1 1 1], 'bgcolor', 'k', 'alpha', 0.1, 'linecolor', 'w', 'LineWidth', 1);
    end
end
 
for rgc = 1:length(mosaic_contours)
    plot_polygon_struct(mosaic_contours{rgc}, 'facecolor', [0 0 0], 'bgcolor', 'w', 'alpha', 0.3, 'linecolor', 'w', 'LineWidth', 1);
end

colormap gray;
