


load '/Volumes/dusom_fieldlab/All_Staff/lab/Development/Natural-Movies/mat-files/bees_mr.mat';

implay(mov, 30)

imagesc(mov(:,:,1))
colormap gray
axis square



% load white noise data
datapath_wn = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-05-13-1/data000/data000';
datarun_wn = load_data(datapath_wn);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_ei(datarun_wn, 'all');
datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
datarun_wn = get_sta_summaries(datarun_wn, 'all');


cell_id = 5266;
temp_index = get_cell_indices(datarun_wn, cell_id);
plot_rf(datarun_wn,cell_id)

temp_fit = datarun_wn.vision.sta_fits{temp_index};
temp_com = temp_fit.mean * datarun_wn.stimulus.stixel_width;

temp_com_absolute = temp_com + [16, 20];


num_frames = floor(7*60.35/2);
background = zeros(800, 600, 'uint8');
for fm = 1:num_frames
    
    temp_img = background;
    temp_img(101:700,:) = squeeze(mov(:,:,1));
    
    imagesc(temp_img)
    colormap gray
    hold on
    [X,Y] = drawEllipse([temp_com_absolute temp_fit.sd*datarun_wn.stimulus.stixel_width temp_fit.angle]);
    coord_tform = coordinate_transform(datarun_wn,'sta');
    [X, Y] = tformfwd(coord_tform, X, Y);
    plot(X,Y,'Color','y', 'LineWidth', 4);
    
end




  [X,Y] = drawEllipse([ctr rad the_fit.angle]);
        
        % if any are NaN, skip
        if any(isnan([X Y]));continue;end
        
        % transform to desired coordinates
        [X, Y] = tformfwd(coord_tform, X, Y);

        % plot the points and fill
        
        if ~strcmpi(params.fill_color, 'none')
            h = fill(X,Y,params.fill_color);
            
            set(h,'facealpha',params.face_alpha)
        end
       plot(X,Y,'Color','k', 'LineWidth', 1);
        hold on
