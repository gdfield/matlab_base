function rf_stat_struct = get_rf_stats(datarun, cell_type, varargin)
% get_rf_stats     This function calculates the mean, ste, and std of the rf radii 
%                       for the cell types listed in cell_type
%
% usage:  [rf_mean_radius, rf_ste_radius, rf_std_radius] = get_rf_stats(datarun, cell_type, varargin)
%
% arguments:     datarun - datarun structure
%               cell_type - names of cell types to use
%            
%
% outputs:      rf_stat_struct.
%               rf_mean_radius - mean radius for cell_type, if cell_type
%                               lists multiple types, rf_mean_radius is a vector
%               rf_ste_radius - ste radius for cell_type, if cell_type
%                               lists multiple types, rf_mean_radius is a vector
%               rf_std_radius - std radius for cell_type, if cell_type
%                               lists multiple types, rf_mean_radius is a vector
%
% optional params, their default values, and what they specify:
%
% fits_to_use           vision          'vision' or 'obivus' fits can be specified
% units                 pixels          'pixels of 'microns' can be specified
% microns_per_pixel       5            For results to be trusted absolutely, this value needs to be
%                                        set according to the datarun
%
% 2017-05 GDF
%

p = inputParser;

% specify list of optional parameters
p.addParameter('fits_to_use', 'vision');
p.addParameter('units', 'pixels');
p.addParameter('microns_per_pixel', 5);
p.addParameter('stixel_thresh', 4.5, @isnumeric)
p.addParameter('area_method', 'hull')
p.parse(varargin{:});

num_types = length(cell_type);
rf_stat_struct.rf_mean_radius = zeros(1, num_types);
rf_stat_struct.rf_ste_radius = zeros(1, num_types);
rf_stat_struct.rf_std_radius = zeros(1, num_types);
rf_stat_struct.mean_rf_area = zeros(1, num_types);
rf_stat_struct.ste_rf_area = zeros(1, num_types);
rf_stat_struct.std_rf_area = zeros(1, num_types);

for rgc_type = 1:num_types

    % fits
    rf_radii = get_rf_fit_radius(datarun, cell_type{rgc_type}, 'fits_to_use', p.Results.fits_to_use, 'units', p.Results.units, 'microns_per_pixel', p.Results.microns_per_pixel) * datarun.stimulus.stixel_height;
    rf_stat_struct.rf_mean_radius(rgc_type) = mean(rf_radii);
    rf_stat_struct.rf_ste_radius(rgc_type) = std(rf_radii) ./ sqrt(length(rf_radii)-1);
    rf_stat_struct.rf_std_radius(rgc_type) = std(rf_radii);
    
    % nonparametric from marks
    rf_areas = get_rf_areas_from_marks(datarun, cell_type{rgc_type}, 'area_method', p.Results.area_method, 'rfs_to_use', 'summaries_filt', 'stixel_thresh', p.Results.stixel_thresh);
    rf_stat_struct.mean_rf_area(rgc_type) = mean(rf_areas);
    rf_stat_struct.ste_rf_area(rgc_type) = std(rf_areas) ./ sqrt(length(rf_areas) - 1);
    rf_stat_struct.std_rf_area(rgc_type) = std(rf_areas);

end

