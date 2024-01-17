function tc_stat_struct = get_tc_stats(datarun, cell_type, varargin)
% get_tc_stats      This function extracts time-to-peak, time-to-zero, and biphasic index 
%                   from STA time course
%
% usage:  tc_stat_struct = get_tc_stats(datarun, cell_type, varargin)
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
% units                 pixels          'pixels of 'microns' can be specified
% ms_per_frame       2/120            For results to be trusted absolutely, this value needs to be
%                                        set according to the datarun
%
% 2017-05 GDF
%

p = inputParser;

% specify list of optional parameters
p.addParameter('units', 'ms');
p.addParameter('ms_per_frame', (2/120)*1000);
p.parse(varargin{:});


num_types = length(cell_type);

% time to peak
tc_stat_struct.mean_peak = zeros(1, num_types);
tc_stat_struct.ste_peak = zeros(1, num_types);
tc_stat_struct.std_peak = zeros(1, num_types);
% time to zero
tc_stat_struct.mean_ttz = zeros(1, num_types);
tc_stat_struct.ste_ttz = zeros(1, num_types);
tc_stat_struct.std_ttz = zeros(1, num_types);
% degree of transience
tc_stat_struct.mean_dot = zeros(1, num_types);
tc_stat_struct.ste_dot = zeros(1, num_types);
tc_stat_struct.std_dot = zeros(1, num_types);
% time to trough
tc_stat_struct.mean_ttt = zeros(1, num_types);
tc_stat_struct.ste_ttt = zeros(1, num_types);
tc_stat_struct.std_ttt = zeros(1, num_types);


for rgc_type = 1:num_types

    temp_indices = get_cell_indices(datarun, cell_type{rgc_type});
    num_rgcs = length(temp_indices);
    
    ttp = zeros(1, num_rgcs);
    tc_dot = zeros(1, num_rgcs);
    ttz = zeros(1, num_rgcs);
    ttt = zeros(1, num_rgcs);
    for rgc = 1:num_rgcs
        temp_tc = datarun.stas.time_courses{temp_indices(rgc)};
        
        % get time to peak
        [~,temp_peak_frame] = ext(temp_tc);
        temp_ttp = temp_peak_frame * p.Results.ms_per_frame;
        ttp(rgc) = temp_ttp;
        
        % get time to zero
        temp_fit = datarun.matlab.sta_fits{temp_indices(rgc)};
        fit_params(1) = temp_fit.scale_one;
        fit_params(2) = temp_fit.scale_two;
        fit_params(3) = temp_fit.tau_one;
        fit_params(4) = temp_fit.tau_two;
        fit_params(5) = temp_fit.n_one_filters;
        fit_params(6) = temp_fit.n_two_filters;
        
        time_range = 4:0.01:10;
        tc_fit = fit_time_course_function(fit_params, time_range);
        [~, temp_ttz] = min(abs(tc_fit));
        ttz(rgc) = ((temp_ttz * 0.01)+4) * (p.Results.ms_per_frame);

        % get DOT
        tc_dot(rgc) = 1-(abs(sum(temp_tc)) ./ sum(abs(temp_tc)));       

        % get time to trough
        time_range = 4:0.01:18;
        tc_fit = fit_time_course_function(fit_params, time_range);

        if ext(temp_tc) < 0
            [~, ttt_frame] = max(tc_fit);
            ttt(rgc) = ((ttt_frame * 0.01)+4) * (2000/120);
        else
            [~, ttt_frame] = min(tc_fit);
            ttt(rgc) = ((ttt_frame * 0.01)+4) * (2000/120);
        end
        
    end
    
    tc_stat_struct.mean_peak(rgc_type) = mean(ttp);
    tc_stat_struct.ste_peak(rgc_type) = std(ttp) ./ (sqrt(num_rgcs) -1);
    tc_stat_struct.std_peak(rgc_type) = std(ttp);  
 
    tc_stat_struct.mean_ttz(rgc_type) = mean(ttz);
    tc_stat_struct.ste_ttz(rgc_type) = std(ttz) ./ (sqrt(num_rgcs) -1);
    tc_stat_struct.std_ttz(rgc_type) = std(ttz);      
    
    tc_stat_struct.mean_dot(rgc_type) = mean(tc_dot);
    tc_stat_struct.ste_dot(rgc_type) = std(tc_dot) ./ (sqrt(num_rgcs) -1);
    tc_stat_struct.std_dot(rgc_type) = std(tc_dot);
    
    tc_stat_struct.mean_ttt(rgc_type) = mean(ttt);
    tc_stat_struct.ste_ttt(rgc_type) = std(ttt) ./ (sqrt(num_rgcs) -1);
    tc_stat_struct.std_ttt(rgc_type) = std(ttt);
    
     
end








