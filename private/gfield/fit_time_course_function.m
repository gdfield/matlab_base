function tc = fit_time_course_function(params, t_points)

%t_points = (1:1:params(7))-1;
t_points;
t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
t_filter_two = params(2) .* (t_points ./ params(4)).^params(6) .* exp(-params(6)*((t_points ./ params(4)) - 1));
tc = t_filter_one + t_filter_two;
%tc = tc(length(t_points):-1:1);