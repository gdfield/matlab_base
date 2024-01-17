function [time_courses] = get_time_courses_matrix(datarun, cell_spec, varargin)
%
% get_time_courses     get time course vectors of several cells and puts
%                      them into a matrix. Note, only works for BW STAs.
%
% usage:  get_time_courses(datarun, cell_ids)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_ids for options)
%
% optional
%       norm_flag       false       if true, set norm of each tc to 1.
%       marks_params    []          see significant_stixels
%
% output:   time_courses: a TxC matrix with T time points and C cells
%
% 2012-08  sneha
% 2017-06   fixed a whole bunch of problems: GDF.

p = inputParser;
p.addParameter('norm_flag', false, @islogical);
p.addParameter('marks_params', [])

p.parse(varargin{:});


cell_indices = get_cell_indices(datarun, cell_spec);
time_courses = zeros(length(cell_indices), datarun.stas.depth);


for rgc = 1:length(cell_indices)
    
    % check to see if time course has already been calculated
    if ~isempty(datarun.stas.time_courses{cell_indices(rgc)})
        time_courses(rgc, :) = datarun.stas.time_courses{cell_indices(rgc)};

    % if not, calculate it
    else
        % check to see if significant stixels have been calculated.
        if isempty(datarun.stas.marks{cell_indices(rgc)})
            
            % get sta
            sta = datarun.stas.stas{cell_index(rgc)};
            % get significant stixels
            sig_stixels = significant_stixels(sta, 'marks_params', marks_params);
        else
            sig_stixels = datarun.stas.marks{cell_indices(rgc)};
        end
        
        tc = time_course_from_sta(sta, sig_stixels);
        time_courses(rgc,:) = tc;
    end
   
     if p.Results.norm_flag
        time_courses(rgc,:) = time_courses(rgc,:) ./ norm(time_courses(rgc,:));
     end
    
end



