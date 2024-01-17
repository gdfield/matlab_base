function opt_model = ln_fit(fitspikes, fitmovie, center, varargin)
%% DESCRIPTION
% NJB 2016-06-03
%
% This takes in the stimulus, the spikes, and the location of the cell to
% fit a LN with a fixed spatial filter and a sigmoidal nonlinearity. 
%
% INPUTS
%
%   REQUIRED
%
%   fitspikes: the spike times of the neuron (in seconds)
%       THESE MUST ALREADY BE ALIGNED TO THE TRIGGERS.
%       Use STA_Test to make sure you can get an STA. Then you will know
%       that they are properly aligned. You can use
%       align_spikes_triggers.m to help with this. 

%   fitmovie: the movie frame by frame. You should
%       have a frame for every 1/monitor_refresh seconds, so if the interval was two, your
%       fitmovie should have 2 of each frame
%       OR the xml specification, like RGB-8-1-0.48-11111-32x32

%   center_coord: the center of the RF. This can be the output from
%       STA_Test
%
%   OPTIONAL
%
%   WN_STA: optional, To do fixedSP_rk1, you need to input the STA in the same
%   dimensions as the fitting stimulus
%
%   monitor_refresh: usually should be 120Hz. This is NOT the interval!
%   Just the monitor speed! This SHOULD be calibrated from the triggers
%
%
% PATHS NEEDED
% Vision.jar, such as javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar')
% The lab codebase, addpath(genpath('Repo location /matlab/code/lab'))
% The glm code folder, addpath(genpath('Repo location /matlab/code/projects/glm))

% MODEL DETAILS
% This code works by first fitting a GLM, and then iteratively fitting a
% temporal filter and a sigmoidal nonlinearity. The spatial filter will be
% fixed as the output of the original GLM fit, so that can be either a
% fixed filter from the STA, or a fit filter. You can set this by changing
% glm_parameters. 
%
% RUN NOTES / TROUBLESHOOTING
% I highly recommend that you include a monitor refresh that is calibrated
% from the triggers. Estimating does not work well for modelling. 
%
% The only difference in inputs is that this DOES NOT take in neighbor
% spikes. The ln model does not include coupling!!!
%
% CHECK THAT YOU CAN GET AN STA USING STA_Test and your fitspikes and
% fitmovie. If you cannot, you need to properly align the spikes to the
% triggers !! Use align_spikes_triggers.m. 
%
% BLOCKED stimuli: If you have interleaved data, rather than one long run,
% use concat_movie and concat_spikes, check that you can
% get an STA using STA_Test, and THEN use this code. 
%
% This also only works for greyscale movies! It is not set up for color!
%
% Check the model architecture in glm_parameters.m
%
% The iterative procedure is NOT OPTIMAL! Should either do more iterations
% or fit together or some criteria for deciding iteration number

%% FUNCTION TIME

% Checks the structure. Obviously this would be more elegant if I changed
% GLMPars for you but it's too annoying so deal.
GLMPars = glm_parameters;
if GLMPars.PostSpikeFilter || GLMPars.CouplingFilters
        error('Set GLMPars.PostSpikeFilter and GLMPars.CouplingFilters in glm_parameters.m in the glm_util folder to false. Just do it.')
end

% Get optional inputs
p = inputParser;
p.addParameter('WN_STA', 0)
p.addParameter('monitor_refresh', 120)
p.addParameter('iterations', 3);
p.addParameter('GLM_params', 'glm_parameters')
p.addParameter('exp', 0);
p.parse(varargin{:});
WN_STA = p.Results.WN_STA;
monitor_refresh = p.Results.monitor_refresh;
iterations = p.Results.iterations;
GLM_params = p.Results.GLM_params;
use_exp = p.Results.exp;
clear p

% More warnings
if monitor_refresh == 120 || monitor_refresh == 119.5
    warning('Calibrating the monitor_refresh using the triggers is a good idea!!')
end

% Fit GLM to get inits and stim details and all that jazz
disp('fitting original GLM')
fittedGLM     = glm_fit(fitspikes, fitmovie,center, 'monitor_refresh', monitor_refresh, 'WN_STA', WN_STA, 'GLM_params',GLM_params );

% fit NL and time iteratively
disp('fitting NL and time')
if iterations == 0
    [model, time] = fit_NL_and_time_noniterative(fittedGLM, fitspikes, fitmovie);
else
    [model, time, iterations] = fit_NL_and_time(fittedGLM, fitspikes, fitmovie, 'iterations', iterations);
    opt_model.iterations = iterations;
end
% structure that shit
opt_model.orig_glm = fittedGLM;
opt_model.model = model;
try 
    opt_model.new_NL.firing_rate = predict(model, [-2:0.5:2]');
catch
    opt_model.new_NL.firing_rate = feval(model, [-2:0.5:2]');
end
opt_model.new_NL.gen_signal = [-2:0.5:2]';
opt_model.new_time = time;


end