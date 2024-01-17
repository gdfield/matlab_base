function [time_to_zero, time_to_peak] = get_time_to_zero(datarun, cell_spec, varargin)

p = inputParser;

% specify list of optional parameters
p.addParameter('units', 'ms');
p.addParameter('ms_per_frame', 2/120);
p.addParameter('time_range', 3:0.01:12);

p.parse(varargin{:});

cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);
time_to_zero = zeros(1, num_rgcs);
time_to_peak = zeros(1, num_rgcs);

for rgc = 1:num_rgcs

    % get time to zero
    temp_fit = datarun.matlab.sta_fits{cell_indices(rgc)};
    fit_params(1) = temp_fit.scale_one;
    fit_params(2) = temp_fit.scale_two;
    fit_params(3) = temp_fit.tau_one;
    fit_params(4) = temp_fit.tau_two;
    fit_params(5) = temp_fit.n_one_filters;
    fit_params(6) = temp_fit.n_two_filters;
    tc_fit = fit_time_course_function(fit_params, p.Results.time_range);
    
     plot(tc_fit)
     pause

    [~, temp_ttz] = min(abs(tc_fit)); 
    time_to_zero(rgc) = ((temp_ttz * 0.01) + p.Results.time_range(1)) * p.Results.ms_per_frame;
    
%     tmp_tc = datarun.stas.time_courses{cell_indices(rgc)};
%     tmp_val = min(tmp_tc ./ ext(tmp_tc)); 
%     trough_val(rgc) = tmp_val;

    [~, tmp_peak] = max(abs(tc_fit));
    time_to_peak(rgc) = ((tmp_peak * 0.01) + p.Results.time_range(1)) * p.Results.ms_per_frame;
end