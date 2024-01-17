function error = timecourse_error_function_NI(timecourse, firing_rate, gen_signal_spatial, mu)
n_bins = 50;

gen_signal = conv(gen_signal_spatial, timecourse, 'full')+mu;
gen_signal = gen_signal(1:(length(firing_rate)));

% Sort the generator signal into percentiles and calculate both the bin
% center and the average firing rate for that bit.
[sorted_gen_signal, index] = sort(gen_signal); % sort the generator signal by valueuse_exp
bin_size = floor(length(gen_signal)/n_bins); % calculate the number of points in each bin
idx = 1:bin_size;
bin = zeros(n_bins,1); bin_FR = zeros(n_bins,1); % initialize those matrices
for i_bin = 1:n_bins
    bin(i_bin) = mean(sorted_gen_signal(idx)); % the center of the bin (in units of the generator signal)
    bin_FR(i_bin) = mean(firing_rate(index(idx))); % the average firing rate for that bin
    idx = idx + bin_size; % move on to the next bin
end

% Fit the NL from the generator signal
%{
if use_exp
    model = fitnlm(bin, bin_FR, 'y~b3+b2*exp(b1*x1)', [0 1 1]);
else
%}
model = fit(bin,bin_FR, 'linearinterp');
%model = fitnlm(bin, bin_FR, 'y~b4+b1/(b2+exp(b3*x1))', NL_init); % Use matlab's fun model fitting to fit the sigmoidal NL

try
    predicted_firing = predict(model, gen_signal');
catch
    predicted_firing = feval(model, gen_signal);
end
%predicted_firing = model_params(1)./(model_params(2)+exp(model_params(3)*gen_signal));
error = sum((predicted_firing' - firing_rate).^2)+10^4*norm(diff(timecourse)); % MSE
end

% fittedGLM.linearfilters.Stimulus.time_rk1'