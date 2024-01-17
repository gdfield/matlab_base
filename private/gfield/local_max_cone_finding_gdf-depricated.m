%data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
%datarun = load_data(data_path);
%datarun = load_sta(datarun,'load_sta',[]);
%new_datarun = datarun;

% data paths
data_path = '/snle/lab/Experiments/Array/Analysis/2011-12-13-2/streamed/data004-0/data004-0';

datarun = load_data(data_path);
datarun = load_params(datarun);

cell_types = {3,4,5};

datarun = load_sta(datarun,'load_sta',cell_types);
datarun = get_sta_summaries(datarun, cell_types);


clear spat_sens_params

% how to combine RGB values of the RF
spat_sens_params.strength = 'vector length';
%spat_sens_params.strength = {'inner or',...
%   [0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};

% how to filter the RF before looking for significant stixels
spat_sens_params.filter =  struct('type','given','filt',...
           make_gaussian('center_radius',0.5,'x_size',5,'y_size',5,'center',[3 3]));

% how to find significant stixels
spat_sens_params.selection_params = struct('type','thresh','thresh',5);
%spat_sens_params.selection_params = struct('type','lasso');
%spat_sens_params.selection_params = struct('type','max','thresh',5);

% how to combine stixels between cells
spat_sens_params.combine_stixels = 'sum';
%cone_loc_params.combine_stixels = 'max';

% online readout of what's going on
spat_sens_params.verbose = false;
spat_sens_params.fig_single_cell = [];
spat_sens_params.foa_spat_sens = 105;


% accumulate the spatial sensitivity across all STAs
[spatial_sensitivity,all_sig_stixels,spatial_cell_ids] = compute_spatial_sensitivity(datarun, cell_types, spat_sens_params);


figure(15)
imagesc(spatial_sensitivity)
axis square
colormap gray
print(15,'~/Desktop/sensitivity2.pdf', '-dpdf')


cone_centers = find_local_maxima(spatial_sensitivity);

% find indices to get cone locations
num_cones = length(find(cone_centers == 1));
[x_indices, y_indices] = find(cone_centers == 1);
centers = zeros(num_cones,2);
centers = [y_indices, x_indices];

%%


% GET PRECISE LOCATION AND COLOR PROFILE OF EACH CONE
% use individual STAs to refine the location of cones and extract their color profiles

% fit BW version of initial cone locations, no penalty

which_cones = 'all';
remap.fcn = @(x)([x(:,1)./x(:,2) x(:,3)./x(:,2)]);
remap.x_caption = 'red/green';
remap.y_caption = 'blue/green';
cones_labeled = bwlabel(spatial_sensitivity,8);

% cone kernel parameters
% this shape is fitted to each cone
cone_kernel_params.type = 'dog';
cone_kernel_params.center_radius = 0.75;
cone_kernel_params.surround_scale = 0;
cone_kernel_params.surround_radius = 1;


% online readout of what's going on
cone_rf_params.regress = 1;
cone_rf_params.combine = 'sum';
cone_rf_params.single_cone_figure = [];%105;
cone_rf_params.cone_weights_figure = [];%106;
cone_rf_params.verbose = true;


initial_cone_centers = cone_centers;

switch which_cones
    case 'all'  % use all cones
        cone_analysis_roi = [];

    case 'roi'  % create a ROI if it doesn't already exist
        if ~exist('cone_analysis_roi','var')
            figure;clf;cone_analysis_roi = roipoly(spatial_sensitivity);
        end
end

% use roi
cone_rf_params.roi = cone_analysis_roi;

% pass in cone remapping function (used only if points are plotted)
cone_rf_params.cone_remap = remap;

% compute cone locations and colors
[cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] =...
    summarize_cone_rfs(datarun, spatial_cell_ids, initial_cone_centers,...
    cones_labeled, all_sig_stixels, cone_kernel_params, cone_rf_params);

% get list of fitted cone centers
cone_centers = zeros(length(cone_ideal_shapes),2);
for cc = 1:length(cone_ideal_shapes)
    cone_centers(cc,1:2) = cone_ideal_shapes{cc}.center;
end


%% plot cone location
y = 300;
x = 300;
rgb = [.5 .5 .5];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(10);clf;image(plot_mat);axis image; hold on
plot(centers(:,1), centers(:,2), 'k.', 'MarkerSize', 15)

datarun.cones.centers = centers;
datarun.cones.types(length(datarun.cones.centers(:,1))) = 'U';

%% compare local max to bayesian cone finding

alt_datarun = datarun;
alt_datarun = load_cones(alt_datarun);

plot_cone_mosaic(alt_datarun, 'fig_or_axes', 10, 'clear', false', 'cone_colors', [1 1 1], 'cone_size', 6)

%% inspect individual cells

cell_indices = get_cell_indices(datarun, cell_types);
window_size = 20;

for cc = 1:length(cell_indices)
    plot_rf(datarun, datarun.cell_ids(cell_indices(cc)), 'foa', 1)
    
    plot_cone_mosaic(datarun, 'fig_or_axes', 1, 'clear', false', 'cone_colors', [0 1 0], 'cone_size', 16)
    plot_cone_mosaic(alt_datarun, 'fig_or_axes', 1, 'clear', false', 'cone_colors', [1 0 1], 'cone_size', 6)
    
    temp_rf_com = datarun.stas.rf_coms{cell_indices(cc)};
    x_bg = temp_rf_com(1) - window_size;
    x_end = temp_rf_com(1) + window_size;
    y_bg = temp_rf_com(2) - window_size;
    y_end = temp_rf_com(2) + window_size;
    
    axis([x_bg x_end y_bg y_end])
    pause
    
end

    
    
    
    
    
    
    
    

