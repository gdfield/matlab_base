
% 3D view
num_rgcs = 36;
space_factor = 8;
overlap_factor = 2;
rf_sd = space_factor * sqrt(3)/2 / 2 * overlap_factor; % space factor * hex mosaic spacing div 2
buf_factor = 3*rf_sd;

locs = hexagonal_mosaic(num_rgcs, 'jitter_method', 'Gaussian', 'spacing', space_factor);
% zero pad around locs
locs = locs + repmat(buf_factor, size(locs));

x_dim = ceil(max(locs(:)) + buf_factor);
y_dim = x_dim;


for rgc = 1:num_rgcs

    rf{rgc} = make_Gaussian_two_d('center_point_x',locs(rgc,1), 'center_point_y', locs(rgc,2),...
                    'sd_x', rf_sd, 'sd_y', rf_sd, 'x_dim', x_dim, 'y_dim', y_dim);     
end

          

sum_surf = zeros(size(rf{1}));
for rgc = 1:num_rgcs
    sum_surf = sum_surf + rf{rgc};
end

surf(rf{10})
view(-38, 68)

surf(sum_surf)
view(-38, 68)

%% 2D view




