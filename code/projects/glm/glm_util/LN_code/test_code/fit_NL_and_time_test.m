function [model, t_init, increments] = fit_NL_and_time_test(fittedGLM, fitspikes, fitmovie, varargin)
% This takes in a fittedGLM structure and makes it better!!!!
% NJB 2016-06-03
% JK but it does iteratively fit a sigmoidal nonlinearity and a timecourse.
% Boom. 

% INPUTS
% 
% REQUIRED: 
%   fittedGLM: a structure that is the output of glm_fit
%   fitspikes and fitmovie: the same ones used to make the fittedGLM

% OPTIONAL:
%   iterations: I found that this converges fast, and 2 is more than
%   enough, but feel free to iterate to your hearts desire. 
%   NL_init: I set these based on some badass parasols, but they might not
%   be ideal for all situations. 
%   n_bins_generator_signal: The number of bins to split the generator
%   signal into to calculate the static nonlinearity. The generator signal
%   will be split into this many bins with an equal number of data points in
%   each bin.
%
%
% NOTES: The "MU" or TONIC FIRING rate gets absorbed here into the
% nonlinearity. Don't worry about it brah!

p = inputParser;
p.addParameter('iterations', 2)
p.addParameter('NL_init', [20 0.5 -2 2]);
p.addParameter('n_bins_generator_signal', 100);
p.parse(varargin{:});
NL_init = p.Results.NL_init;
iterations = p.Results.iterations;
n_bins = p.Results.n_bins_generator_signal;
clear p
tic

% Initiate the temporal filter with the output of the GLM
t_init = fittedGLM.linearfilters.Stimulus.time_rk1;
increments.NL.start = NL_init;
increments.time.start = t_init;

% Iterate
for iter = 1:iterations
    disp(['iteration ' num2str(iter) ': fitting NL'])
    
    % First calculate the generator signal and firing rate for all the binz
    gen_signal_spatial = glm_gen_signal_spatial(fittedGLM, fitmovie);
    gen_signal = conv(gen_signal_spatial, t_init, 'full'); % MU/TONIC FIRING gets absorbed in the nonlinearity, so all we need is the spatial component
    home_spbins  = ceil(fitspikes / (fittedGLM.bins_per_frame*fittedGLM.t_bin));
    home_spbins = home_spbins(home_spbins < length(gen_signal));
    firing = zeros(size(gen_signal));
    for i=1:length(home_spbins)
        firing(home_spbins(i)) = firing(home_spbins(i))+1;
    end
    window = 10; % I don't think smoothing is necessary
    firing_rate = conv(firing, gausswin(window), 'same')/(sum(gausswin(window))*(fittedGLM.bins_per_frame*fittedGLM.t_bin)); % go from spikes to firing rates
    
    % Sort the generator signal into percentiles and calculate both the bin
    % center and the average firing rate for that bit.
    [sorted_gen_signal, index] = sort(gen_signal); % sort the generator signal by value
    bin_size = floor(length(gen_signal)/n_bins); % calculate the number of points in each bin 
    idx = 1:bin_size;
    bin = zeros(n_bins,1); bin_FR = zeros(n_bins,1); % initialize those matrices
    for i_bin = 1:n_bins
        bin(i_bin) = mean(sorted_gen_signal(idx)); % the center of the bin (in units of the generator signal)
        bin_FR(i_bin) = mean(firing_rate(index(idx))); % the average firing rate for that bin
        idx = idx + bin_size; % move on to the next bin
    end
    
    % Fit the NL from the generator signal
    model = fitnlm(bin, bin_FR, 'y~b4+b1/(b2+exp(b3*x1))', NL_init); % Use matlab's fun model fitting to fit the sigmoidal NL
    NL_init = model.Coefficients.Estimate; % set the initial values for the next iteration
    
    % Fit the time course
    disp(['iteration ' num2str(iter) ': fitting timecourse'])
    new_time = fminunc(@(time_course)timecourse_error_function(time_course, firing_rate, model, gen_signal_spatial, 0), t_init);
    t_init = new_time; % set the initial values for the next iteration
    clear new_time;
    
    % Save the values at each iteration to check convergence post 
    increments.NL.iter{iter} = NL_init;
    increments.time.iter{iter} = t_init;
    toc
end
end