

datapath = '/Users/gfield/Analysis/rat/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/data000-map/data000-map';


%%
cell_types = {'ON type1'};
datarun = load_data(datapath);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', cell_types);
temp_indices = get_cell_indices(datarun, cell_types);
num_rgcs = length(temp_indices);
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, cell_types, 'marks_params', marks_params, 'sta', 'stas');


% KEY FUNCTION HERE
datarun = tc_center_surround_fit(datarun, cell_types);


% This loop plots comparisons between data and fits
for rgc = 1:num_rgcs
    
    tmp_fit_params = datarun.matlab.sta_fits{temp_indices(rgc)};
    fit_params = extract_fit_params_from_fit(tmp_fit_params); 
    temp_tc = datarun.stas.time_courses{temp_indices(rgc)};

    temp_fit_sta = sta_fit_function(fit_params);
    
    % get the significant stixels
    temp_sig_stixels = significant_stixels(temp_fit_sta, 'time', 'max', 'thresh', 4.5);
    temp_fit_time_course = time_course_from_sta(temp_fit_sta, temp_sig_stixels);

    % Assay fit quality
    figure(1); clf;

    % compare time courses
    subplot(2,2,1)
    plot(temp_fit_time_course ./ norm(temp_fit_time_course), 'k')
    hold on
    temp_tcs = datarun.stas.time_courses{temp_indices(rgc)};
    plot(temp_tc ./ norm(temp_tc), 'r')
    %legend('fit', 'data', 'Location', 'northwest')
    title('time course')
    xlabel('frames')
    hold off

    % plot fit RF
    subplot(2,2,2)
    reshape_sta = reshape(squeeze(temp_fit_sta), [], size(temp_fit_sta,4));
    fit_rf = reshape_sta * temp_fit_time_course;
    fit_rf = reshape(fit_rf, 40, 80);
    fit_rf = repmat(fit_rf, [1,1,3]);
    image(norm_image(fit_rf))
    title('fit RF')

    subplot(2,2,3)
    temp_rf = datarun.stas.rfs{temp_indices(rgc)};
    image(norm_image(temp_rf));
    title('RF')

    subplot(2,2,4)
    res_rf = (temp_rf ./ max(temp_rf(:))) - (fit_rf ./ max(fit_rf(:)));
    image(norm_image(res_rf))
    title('residual RF')
    drawnow
    pause
end
