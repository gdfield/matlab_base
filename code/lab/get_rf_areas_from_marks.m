function rf_areas = get_rf_areas_from_marks(datarun, cell_spec, varargin)
% get_rf_areas_from_marks     This function calculates the mean, ste, and std of the rf area 
%                           from the marks
%
% usage:  rf_areas = get_rf_areas_from_marks(datarun, cell_spec, varargin)
%
% arguments:     datarun - datarun structure
%               cell_spec - cell_spec, see get_cell_indices
%            
%
% outputs:     rf_areas - area of all rfs for cells in cell_spec
%
%
% optional params, their default values, and what they specify:
%
% units                 stixels          'pixels of 'microns' can be specified
% microns_per_pixel       25            For results to be trusted absolutely, this value needs to be
%                                        set according to the datarun
%
% 2017-05 GDF


p = inputParser;

% specify list of optional parameters
p.addParameter('micron_sq_per_stixel', 25, @isnumeric);
p.addParameter('area_method', 'hull');
p.addParameter('rfs_to_use', [])
p.addParameter('stixel_thresh', 4.5, @isnumeric)
p.parse(varargin{:});

temp_indices = get_cell_indices(datarun, cell_spec);
rf_areas = zeros(1, length(temp_indices));

for rgc = 1:length(temp_indices)
    
    switch p.Results.area_method
        case 'sig_stixels'
            
            if isempty(p.Results.rfs_to_use)
                rf_areas(rgc) = length(find(datarun.stas.marks{temp_indices(rgc)})) * p.Results.micron_sq_per_stixel;
    
            elseif strcmp(p.Results.rfs_to_use, 'summaries_filt')
                temp_rf = datarun.stas.summaries_filt{temp_indices(rgc)};
                stix = significant_stixels(temp_rf, 'select', 'thresh', 'thresh', p.Results.stixel_thresh);
                rf_areas(rgc) = length(find(stix));

            end

                
        case 'hull'
            if isempty(p.Results.rfs_to_use)
                [tempY, tempX] = find(datarun.stas.marks{temp_indices(rgc)});
                
            elseif strcmp(p.Results.rfs_to_use, 'summaries_filt')
                temp_rf = datarun.stas.summaries_filt{temp_indices(rgc)};
                stix = significant_stixels(temp_rf, 'select', 'thresh', 'thresh', p.Results.stixel_thresh);
                [tempY, tempX] = find(stix);
            end

            if isempty(tempY)
                rf_areas(rgc) = NaN;
                warning(['RGC ', num2str(datarun.cell_ids(temp_indices(rgc))), ' in ', datarun.names.short_name, ' has no sig stixels'])
            else
                [~, temp_area] = convhull(tempX,tempY);
                rf_areas(rgc) = temp_area .* p.Results.micron_sq_per_stixel;
            end
    end
            
end



