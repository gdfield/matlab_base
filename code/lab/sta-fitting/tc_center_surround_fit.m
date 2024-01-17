function datarun = tc_center_surround_fit(datarun, cell_spec, varargin)
% tc_center_surround_fit     This function performs a fit on STAs in data
%       specified by cell spec. The fit uses a serial approach that first
%       fits the temporal RF in isolation, then the spatial RF center, then
%       the spatial RF surround, then uses this as a set of initial conditions
%       on a full spatiotemporal fit. The function utilizes MATLAB's parralel
%       processing toolbox with a parfor loop. 
%
% Note: the code has only been debuggd for BW WN stimuli.
%
% usage:   datarun = tc_center_surround_fit(datarun, cell_spec, varargin)
%
% arguments:     datarun - datarun structure
%               cell_type - names of cell types to use
%            
%
% outputs:      datatun structure fit fits located in
%               dataun.matlab.sta_fits
%
% optional params, their default values, and what they specify:
%
% verbose               false          report cell that is being fit 

%
% 2017-05 GDF
%

% process variable inputs
p = inputParser;
p.addParameter('verbose', true, @islogical);
p.parse(varargin{:});
verbose = p.Results.verbose;

% parse cell types, number of RGCs. and STA size.
cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);
[height, width, color, num_frames] = size(datarun.stas.stas{cell_indices(1)});

% initialize matrixes for TCs and STAs and put info into these matrices
% Note: this is needed for efficient use of parfor function
tc_matrix = zeros(length(datarun.stas.time_courses{cell_indices(1)}), num_rgcs);
sta_matrix = zeros(numel(datarun.stas.stas{cell_indices(1)}), num_rgcs);
for rgc = 1:num_rgcs
    temp_tc = datarun.stas.time_courses{cell_indices(rgc)};
    tc_matrix(:, rgc) = temp_tc;
    temp_sta = datarun.stas.stas{cell_indices(rgc)};
    sta_matrix(:, rgc) = reshape(temp_sta, [],1);
end

% get the RFs from the data run for plotting.
if verbose
    rf_matrix = zeros(numel(datarun.stas.rfs{cell_indices(1)}), num_rgcs);
    for rgc = 1:num_rgcs
        temp_rf = datarun.stas.rfs{cell_indices(rgc)};
        rf_matrix(:,rgc) = reshape(temp_rf, [], 1);
    end
end
cell_ids = datarun.cell_ids;
        
parfor rgc = 1:num_rgcs

    if verbose
        fprintf('fitting the STA for cell %d... \n', cell_ids(cell_indices(rgc)))
    end
 
    
    % get the tc and sta for the rgc
    %temp_tc = datarun.stas.time_courses{temp_indices(rgc)};
    %temp_sta = datarun.stas.stas{temp_indices(rgc)};
    temp_tc = tc_matrix(:,rgc);
    temp_sta = sta_matrix(:,rgc);
    temp_sta = reshape(temp_sta, [height, width, color, num_frames]);
     
    %%%%%%%%%%%%% fit just the time course %%%%%%%%%%%%%
    [~, final_params] = fit_time_course(temp_tc);

    
    %%%%%%%%%%%%% fit the center %%%%%%%%%%%%%
    fit_instructions = [];
    fit_instructions.fit_sig_stixels_only = false;
    fit_instructions.fit_surround = false;
    fit_instructions.fit_surround_sd_scale = false;
    fit_instructions.fit_surround_amp_scale = false;
    fit_instructions.fit_center_point_x = true;
    fit_instructions.fit_center_point_y = true;
    fit_instructions.fit_center_sd_x = true;
    fit_instructions.fit_center_sd_y = true;
    fit_instructions.fit_center_rotation_angle = true;
    fit_instructions.fit_center_amp_scale = true;
    fit_instructions.fit_color_weight_a = true;
    fit_instructions.fit_color_weight_b = true;
    fit_instructions.fit_color_weight_c = true;
    fit_instructions.fit_scale_one = false;
    fit_instructions.fit_scale_two = false;
    fit_instructions.fit_tau_one = false;
    fit_instructions.fit_tau_two = false;
    fit_instructions.fit_n_one_filters = false; % adjust number of lowpass filters

    % set fit params from previous fit
    fit_instructions.initial_scale_one = final_params(1); %scale_one;
    fit_instructions.initial_scale_two = final_params(2); %scale_two;
    fit_instructions.initial_tau_one = final_params(3); %tau_one;
    fit_instructions.initial_tau_two = final_params(4); %tau_two;
    fit_instructions.initial_n_one_filters = final_params(5); %n_one_filters;
    fit_instructions.initial_n_two_filters = final_params(6); %n_two_filters;
   
    % Perform the center only fit
    temp_fit_params = fit_sta(temp_sta, fit_instructions);    

    
    %%%%%%%%%%%%% fit the surround %%%%%%%%%%%%%
    fit_instructions = [];
    fit_instructions.fit_sig_stixels_only = false;
    fit_instructions.fit_surround = true;
    fit_instructions.fit_surround_sd_scale = true;
    fit_instructions.fit_surround_amp_scale = true;
    % hold center and TC
    fit_instructions.fit_center_point_x = false;
    fit_instructions.fit_center_point_y = false;
    fit_instructions.fit_center_sd_x = true;
    fit_instructions.fit_center_sd_y = true;
    fit_instructions.fit_center_rotation_angle = false;
    fit_instructions.fit_center_amp_scale = true;
    fit_instructions.fit_color_weight_a = false;
    fit_instructions.fit_color_weight_b = false;
    fit_instructions.fit_color_weight_c = false;
    fit_instructions.fit_scale_one = false;
    fit_instructions.fit_scale_two = false;
    fit_instructions.fit_tau_one = false;
    fit_instructions.fit_tau_two = false;
    fit_instructions.fit_n_one_filters = false; % adjust number of lowpass filters

    % set fit params from previous fit
    fit_instructions.initial_scale_one = final_params(1); %scale_one;
    fit_instructions.initial_scale_two = final_params(2); %scale_two;
    fit_instructions.initial_tau_one = final_params(3); %tau_one;
    fit_instructions.initial_tau_two = final_params(4); %tau_two;
    fit_instructions.initial_n_one_filters = final_params(5); %n_one_filters;
    fit_instructions.initial_n_two_filters = final_params(6); %n_two_filters;
    
    fit_instructions.initial_center_point_x = temp_fit_params.center_point_x;
    fit_instructions.initial_center_point_y = temp_fit_params.center_point_y;
    fit_instructions.initial_center_rotation_angle = temp_fit_params.center_rotation_angle;
    fit_instructions.initial_color_weight_a = temp_fit_params.color_weight_a;
    fit_instructions.initial_color_weight_b = temp_fit_params.color_weight_b;
    fit_instructions.initial_color_weight_c = temp_fit_params.color_weight_c;
    fit_instructions.initial_surround_sd_scale = 2;
    fit_instructions.initial_surround_amp_scale = 0.666/sqrt(2);
     
    % perform the surround fit
    temp_fit_params_two = fit_sta(temp_sta, fit_instructions);    

   
    %%%%%%%%%%%%% fit ALL %%%%%%%%%%%%%
    fit_instructions = [];
    fit_instructions.fit_sig_stixels_only = false;
    fit_instructions.fit_surround = true;
    fit_instructions.fit_surround_sd_scale = true;
    fit_instructions.fit_surround_amp_scale = true;
    % hold center and TC
    fit_instructions.fit_center_point_x = true;
    fit_instructions.fit_center_point_y = true;
    fit_instructions.fit_center_sd_x = true;
    fit_instructions.fit_center_sd_y = true;
    fit_instructions.fit_center_rotation_angle = true;
    fit_instructions.fit_center_amp_scale = true;
    fit_instructions.fit_color_weight_a = false;
    fit_instructions.fit_color_weight_b = false;
    fit_instructions.fit_color_weight_c = false;
    fit_instructions.fit_scale_one = true;
    fit_instructions.fit_scale_two = true;
    fit_instructions.fit_tau_one = true;
    fit_instructions.fit_tau_two = true;
    fit_instructions.fit_n_one_filters = true; % adjust number of lowpass filters
    fit_instructions.fit_n_two_filters = true; % adjust number of lowpass filters

    % set fit params from previous fit
    fit_instructions.initial_scale_one = final_params(1); %scale_one;
    fit_instructions.initial_scale_two = final_params(2); %scale_two;
    fit_instructions.initial_tau_one = final_params(3); %tau_one;
    fit_instructions.initial_tau_two = final_params(4); %tau_two;
    fit_instructions.initial_n_one_filters = final_params(5); %n_one_filters;
    fit_instructions.initial_n_two_filters = final_params(6); %n_two_filters;
    
    fit_instructions.initial_center_point_x = temp_fit_params_two.center_point_x;
    fit_instructions.initial_center_point_y = temp_fit_params_two.center_point_y;
    fit_instructions.initial_center_rotation_angle = temp_fit_params_two.center_rotation_angle;
    fit_instructions.initial_color_weight_a = temp_fit_params_two.color_weight_a;
    fit_instructions.initial_color_weight_b = temp_fit_params_two.color_weight_b;
    fit_instructions.initial_color_weight_c = temp_fit_params_two.color_weight_c;
    
    fit_instructions.initial_surround_sd_scale = temp_fit_params_two.surround_sd_scale;
    fit_instructions.initial_surround_amp_scale = temp_fit_params_two.surround_amp_scale;

    % perform the 'free' fit
    temp_fit_params_final = fit_sta(temp_sta, fit_instructions);    

%     if 1
%         fit_params = extract_fit_params_from_fit(temp_fit_params_final); 
%         temp_fit_sta_three = sta_fit_function(fit_params);
%         % get the significant stixels
%         temp_sig_stixels = significant_stixels(temp_fit_sta_three, 'time', 'max', 'thresh', 5);
%         temp_fit_time_course = time_course_from_sta(temp_fit_sta_three, temp_sig_stixels);
% 
%         % Assay fit quality
%         figure(1); clf;
% 
%         % compare time courses
%         subplot(2,2,1)
%         plot(temp_fit_time_course ./ norm(temp_fit_time_course), 'k')
%         hold on
%         plot(temp_tc ./ norm(temp_tc), 'r')
%         legend('fit', 'data', 'Location', 'northwest')
%         title('time course')
%         xlabel('frames')
%         hold off
% 
%         % plot fit RF
%         subplot(2,2,2)
%         reshape_sta = reshape(squeeze(temp_fit_sta_three), [], size(temp_fit_sta_three,4));
%         fit_rf = reshape_sta * temp_fit_time_course;
%         fit_rf = reshape(fit_rf, 40, 80);
%         fit_rf = repmat(fit_rf, [1,1,3]);
%         image(norm_image(fit_rf))
%         title('fit RF')
% 
%         subplot(2,2,3)
%         temp_rf = reshape(rf_matrix(:,rgc), [height, width]);
%         image(norm_image(temp_rf));
%         title('RF')
% 
%         subplot(2,2,4)
%         res_rf = (temp_rf ./ max(temp_rf(:))) - (fit_rf ./ max(fit_rf(:)));
%         image(norm_image(res_rf))
%         title('residual RF')
%         drawnow
%     end
    
    final_fits{rgc} = temp_fit_params_final;
    
end

for rgc = 1:num_rgcs
    datarun.matlab.sta_fits{cell_indices(rgc)} = final_fits{rgc};
end


