function nnnds = get_nnnds(datarun, cell_spec, varargin)
%
% usage:  nnnds = get_nnnds(datarun, cell_spec, varargin)
%
%
% outputs:     nnnds - nnnd between  cells
%
% optional params, their default values, and what they specify:
%
% fits_to_use           matlab          'vision' or 'matlab' fits can be specified
%
% 2018-04 GDF


p = inputParser;

% specify list of optional parameters
p.addParameter('fits_to_use', 'matlab');
p.parse(varargin{:});

fits_to_use = p.Results.fits_to_use;


cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);
all_dists = zeros(num_rgcs);

switch fits_to_use
    
    case 'vision'
        sta_fits = datarun.vision.sta_fits;
        
        for rgc1 = 1:num_rgcs
            for rgc2 = 1:num_rgcs

                fit1 = sta_fits{cell_indices(rgc1)};
                fit2 = sta_fits{cell_indices(rgc2)};

                d = nnnd_gdf(fit1.mean, fit2.mean, fit1.angle, fit2.angle, fit1.sd, fit2.sd);
                all_dists(rgc1, rgc2) = d;
            end
        end

    case 'matlab'
        sta_fits = datarun.matlab.sta_fits;

        for rgc1 = 1:num_rgcs
            for rgc2 = 1:num_rgcs
        
                fit1 = sta_fits{cell_indices(rgc1)};
                fit2 = sta_fits{cell_indices(rgc2)};

                com1 = [fit1.center_point_x, fit1.center_point_y];
                angle1 = fit1.center_rotation_angle;
                sigmas1 = [fit1.center_sd_x, fit1.center_sd_y];

                com2 = [fit2.center_point_x, fit2.center_point_y];
                angle2 = fit2.center_rotation_angle;
                sigmas2 = [fit2.center_sd_x, fit2.center_sd_y];

                d = nnnd_gdf(com1, com2, angle1, angle2, sigmas1, sigmas2);
                all_dists(rgc1, rgc2) = d;
            end
        end
        
    % get the nearest neighbors
end

nnnds = zeros(num_rgcs, 1);
for rgc = 1:num_rgcs
    temp_dists = sort(unique(all_dists(rgc,:)), 'ascend');
    nnnds(rgc) = min(temp_dists(2:end));
end


      
     




    
    
        
        
        
        
        
        
        