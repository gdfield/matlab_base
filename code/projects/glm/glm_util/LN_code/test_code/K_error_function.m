function error = K_error_function(combo_filter, firing_rate, model, fittedGLM, fitmovie, mu)
timecourse = combo_filter(1:30);
space = combo_filter(31:end);
% centerx,centery, sd_center, sd surround, angle
spatial{1} = make_Gaussian_two_d('center_point_x', space(1), 'center_point_y', space(2), 'sd_x', space(3), 'sd_y',  space(3), 'rotation_angle',space(5), 'x_dim', fittedGLM.GLMPars.stimfilter.ROI_length, 'y_dim', fittedGLM.GLMPars.stimfilter.ROI_length);
spatial{2} = make_Gaussian_two_d('center_point_x', space(1), 'center_point_y', space(2), 'sd_x', space(4), 'sd_y',  space(4), 'rotation_angle', space(5), 'x_dim', fittedGLM.GLMPars.stimfilter.ROI_length, 'y_dim', fittedGLM.GLMPars.stimfilter.ROI_length);
space = spatial{1}-spatial{2};
fittedGLM_temp = fittedGLM;
fittedGLM_temp.linearfilters.Stimulus.space_rk1 = space;
gen_signal_spatial = glm_gen_signal_spatial(fittedGLM_temp, fitmovie);
gen_signal = conv(gen_signal_spatial, timecourse, 'full')+mu;
gen_signal = gen_signal(1:(length(firing_rate)));
predicted_firing = predict(model, gen_signal');
error = sum((predicted_firing' - firing_rate).^2);%+10^4*norm(diff(timecourse)); % MSE
end

% fittedGLM.linearfilters.Stimulus.time_rk1'