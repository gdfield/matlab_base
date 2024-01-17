function [GLMT, GLMP] = glm_parameters_84
% Please take a careful look through these parameters! They can make or
% break your model fit.


%% MODEL ARCHITECTURE

%% Spatial Filter Type. Choose FixedSP, Rk1 or Rk2
GLMT.stimfilter_mode = 'fixedSP_rk1_linear'; GLMT.CONVEX = true; % The spatial part is taken from the STA
%GLMT.stimfilter_mode = 'rk1'; GLMT.CONVEX = false; % Both the spatial and temporal parts are fit, rank 1
% GLMT.stimfilter_mode = 'rk2'; GLMT.CONVEX = false; % Two spatial and temporal parts are fit, rank 2

%% Coupling on or off?
GLMT.CouplingFilters = false; % If this is on, you need to include neighborspikes as a parameter to glm_fit and glm_predict

%% Use this to do a quick check that things are right
GLMT.debug      = false; % this just increases the tolerance so the fitting stops more quickly. 
% only set to true if you are testing to make sure everything runs and
% you get a semi-reasonable result. Set to false for your real model fits!

%% Do you want a post spike filter?
GLMT.PostSpikeFilter = true;

%% You probably don't want to change these things
GLMT.cone_model = '8pix_Identity_8pix'; GLMT.cone_sname='p8IDp8';
GLMT.TonicDrive = true; % background firing rate term on or off?
GLMT.StimFilter = true;
GLMT.nullpoint  = 'mean';
GLMT.map_type   = 'mapPRJ';

%% You definitely don't want to change these things
 % these are all for testing and probably don't work. leave set to false. 
GLMT.specialchange = false;
GLMT.Subunits = false;
GLMT.contrast = false;
GLMT.init = false;

%% GLM PARAMETERS
% Detailed parameters of the model, including time bins, filter bases, etc.
% You probably don't want to mess with these much, but if you want to tweak
% the model fits, here ya go.

%% Timing
GLMP.bins_per_frame = 40;
%GLMP.approx_refresh_hz  = 120;
%GLMP.approx_refresh_hz  = 119.5; removing these parameters because why
GLMP.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMP.timenotes_1    = 'True tstim is usually .0083275' ;
GLMP.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMP.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';

%% STA and Spatial Filter Size

GLMP.stimfilter.fixedSP_type = 'WNSTA';
GLMP.stimfilter.ROI_length = 13;  % This is the SIZE of one dimension of the spatial filter in pixels
if ~mod(GLMP.stimfilter.ROI_length, 2)
    error('Please use an ODD ROI length for indexing purposes. JUST DO IT. Edit in this function (glm_parameters)')
end
GLMP.stimfilter.frames = 15;  % This is the number of frames in the temporal filter

%% Probably don't need to worry about the rest of this stuff, although if you are getting weird post spike filter fits,
% you can change the ps basis parameters 

GLMP.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMP.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMP.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMP.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';

% Optimization
GLMP.optimization.tolfun   = 5;
GLMP.optimization.tolx     = 9;
GLMP.optimization.maxiter  = 300; 
GLMP.optimization.note1    = 'tolfun: significant digits of MLE function being optimized';
GLMP.optimization.note2    = 'tolx: significant digits of the input variables of the MLE function';

% Post-spike and Coupling Basis settings
GLMP.spikefilters.ps_note = 'parameters regarding the post-spike filter';
GLMP.spikefilters.cp_note = 'parameters regaring coupling filters';
GLMP.spikefilters.note0 = 'all parameters related to raised sinusoidal humps';
GLMP.spikefilters.note1 = 'basis built by prep_spikefilterbasisGP / create_histbasis as of 2014-05-3';
GLMP.spikefilters.ps.ms  = 100 ;      %% post spike filter time length in millisecs
GLMP.spikefilters.cp.ms  = 100 ;      %% cp spike filter time length in millisecs
%GLMPars.spikefilters.spcng_psf = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
%GLMPars.spikefilters.spcng_cp  = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMP.spikefilters.BiDirect_CP     = false;
GLMP.spikefilters.ps.filternumber = 8;
GLMP.spikefilters.cp.filternumber = 8;
GLMP.spikefilters.ps.spacing      = pi/4;
GLMP.spikefilters.cp.spacing      = pi/2;
GLMP.spikefilters.ps.bstretch     = 0.95;
GLMP.spikefilters.ps.alpha        = 0;
GLMP.spikefilters.cp.bstretch     = 0.05;
GLMP.spikefilters.cp.alpha        = 0;
GLMP.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMP.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out
GLMP.spikefilters.cp.n_couplings = 6  ;
GLMP.others.fitblockchange = false;
end


