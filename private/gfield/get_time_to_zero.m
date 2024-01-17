function [time_to_zero, time_to_peak] = get_time_to_zero(datarun, cell_spec, varargin)

p = inputParser;

% specify list of optional parameters
p.addParameter('ms_per_frame', 2/120);
p.addParameter('step_size', 0.01);
p.addParameter('start_time', 3);
p.addParameter('end_time', 12);

p.parse(varargin{:});

cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);
time_to_zero = zeros(1, num_rgcs);
time_to_peak = zeros(1, num_rgcs);

time_samples = p.Results.start_time:p.Results.step_size:p.Results.end_time;

temp_time = p.Results.start_time:p.Results.step_size:p.Results.end_time;
temp_time = temp_time * p.Results.ms_per_frame;


for rgc = 1:num_rgcs

    % get time to zero
    temp_fit = datarun.matlab.sta_fits{cell_indices(rgc)};
    fit_params(1) = temp_fit.scale_one;
    fit_params(2) = temp_fit.scale_two;
    fit_params(3) = temp_fit.tau_one;
    fit_params(4) = temp_fit.tau_two;
    fit_params(5) = temp_fit.n_one_filters;
    fit_params(6) = temp_fit.n_two_filters;
    tc_fit = fit_time_course_function(fit_params, time_samples);
    

     

    % compute time to zero
    [~, temp_ttz] = min(abs(tc_fit)); 
    time_to_zero(rgc) = ((temp_ttz * p.Results.step_size) + p.Results.start_time) * p.Results.ms_per_frame;

    if time_to_zero(rgc) > 0.3 || time_to_zero(rgc) < 0.1
        time_to_zero(rgc) = NaN;
%         extended_samples = p.Results.start_time:p.Results.step_size:p.Results.end_time* 2;
%         tc_fit = fit_time_course_function(fit_params, extended_samples);
%         [~, temp_ttz] = min(abs(tc_fit)); 
%         time_to_zero(rgc) = ((temp_ttz * p.Results.step_size) + p.Results.start_time) * p.Results.ms_per_frame;
%     end
%     if time_to_zero(rgc) == p.Results.end_time *2
%         time_to_zero(rgc) = NaN;
    end
    
%     plot(temp_time, tc_fit)
%     time_to_zero(rgc)
%     pause
    
    % compute time to trough
%     tmp_tc = datarun.stas.time_courses{cell_indices(rgc)};
%     tmp_val = min(tmp_tc ./ ext(tmp_tc)); 
%     trough_val(rgc) = tmp_val;

    % compute time to peak
    [~, tmp_peak] = max(abs(tc_fit));
    time_to_peak(rgc) = ((tmp_peak * p.Results.step_size) + p.Results.start_time) * p.Results.ms_per_frame;
end