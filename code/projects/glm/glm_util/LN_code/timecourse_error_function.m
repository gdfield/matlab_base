function error = timecourse_error_function(timecourse, firing_rate, model, gen_signal_spatial, mu)
gen_signal = conv(gen_signal_spatial, timecourse, 'full')+mu;
gen_signal = gen_signal(1:(length(firing_rate)));
try
    predicted_firing = predict(model, gen_signal');
catch
    predicted_firing = feval(model, gen_signal);
end
%predicted_firing = model_params(1)./(model_params(2)+exp(model_params(3)*gen_signal));
error = sum((predicted_firing' - firing_rate).^2)+10^4*sum(diff(timecourse).^2); % MSE
end

% fittedGLM.linearfilters.Stimulus.time_rk1'