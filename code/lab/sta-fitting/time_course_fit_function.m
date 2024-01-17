function tc = time_course_fit_function(params)
%
% Usage: tc_fit = time_course_fit_function(params)
%
% Inputs: 
%   params              A vector of parameters with the order ascribed below   
%
% Outputs:
%   tc                  a time course vector
%
% This low level function receives a list of params (the order of which is
% important) and makes a time course based on those
% parameters.
%
%
% PARAM INFORMATION
% each number is an index into the vector 'params' and the name next to the
% number indicates how that value is used.
%
% 1 scale_one
% 2 scale_two
% 3 tau_one
% 4 tau_two
% 5 n-filters-1
% 6 n-filters-2
% 7 frame_number
%
% Author: GDF 
% Data: 2014-05
%


% BODY OF FUNCTION

% time course
%%%%
t_points = (1:1:params(7))-1;
t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
t_filter_two = params(2) .* (t_points ./ params(4)).^params(6) .* exp(-params(6)*((t_points ./ params(4)) - 1));
tc = t_filter_one + t_filter_two;
tc = tc(params(7):-1:1);


