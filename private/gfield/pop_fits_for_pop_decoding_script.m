%% Field Lab multi-neuron model fits -- initial demonstration

addpath(genpath('~/Documents/MATLAB/NIMclass-master'))
addpath(genpath('~/Documents/MATLAB/markSchmidt'))

cd ~/Desktop/
%load white_noise_example_data.mat
%load repeat_data.mat
load NIM_decoding_test.mat

white_noise_spikes = OffBriskTransient_spikes;

% Establish experimental parameters
NC = length(white_noise_spikes);
[NX, NY, NT] = size(white_noise_movie);
% organize stim to be +/- 1
stim = reshape(white_noise_movie,[NX*NY,NT])';
stim(stim < 0.5) = -1;
stim(stim > 0.5) = 1;

% reduce stim length to speed fitting and/or stay in memory limits
num_frames = size(stim,1);
cut_time = (num_frames * refresh) + white_noise_movie_start;
stim = stim(1:num_frames,:);

% % Bin spike train
frac = 8;
Robs = zeros(frac*num_frames, NC); % third version at 8x-resolution (2 ms)
tbins = (0:frac*num_frames)*refresh/frac+white_noise_movie_start;

for rgc = 1:NC
	dum = histc(white_noise_spikes{rgc}(white_noise_spikes{rgc} < cut_time), tbins);
	Robs(:,rgc) = dum(1:end-1);
end


%% Model setup
lat_skip = 0;  % how many time lags to throw out preceding stimulus
num_lags = 20; % number of time lags in temporal kernel: note that response is analyzed here at refresh/2

% prep fitting parameters
stim_par = NIM.create_stim_params( [num_lags NX NY], 'stim_dt', refresh, 'upsampling', frac, 'tent_spacing', frac );
Xstim = NIM.create_time_embedding( stim, stim_par ); % design stim matrix
Rsh = NIM.shift_mat_zpad( Robs, -lat_skip );  % shifted responses
[Ui, Xi] = NIM.generate_XVfolds( num_frames*frac ); % establish 5-fold cross-validation indices


% STAs
stas = Xstim'*Rsh;
figure
for rgc = 1:NC
	subplot(3,4, rgc)
	k = reshape(stas(:,rgc), [num_lags, NX*NY]);
	plot(k)
end

fit_NIM_flag = true; % if true, fit NIM, else only fit GLM
silent = false;
% GLM
tic
clear fit_struct
for rgc = 1:NC
    %fit each RGC
    disp(['%%%% FITTING RGC: ', num2str(rgc), ' %%%'])

    LN_fit = NIM( stim_par, {'lin'}, 1, 'd2xt', 100 );
    LN_fit = LN_fit.fit_filters( Rsh(:,rgc), Xstim, Ui, 'silent', silent );
    %glm0 = glm0.reg_path( Rsh(:,rgc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' );

    % Add spike-history term
    glm_fit = LN_fit.init_spkhist( 12, 'doubling_time', 4, 'negcon');
    glm_fit = glm_fit.fit_filters( Rsh(:,rgc), Xstim, Ui, 'silent', silent );
    %glm_fit.display_model('Xstims', Xstim);
    %drawnow

    fit_struct(rgc).LN_fit = LN_fit;
    fit_struct(rgc).glm_fit = glm_fit;
    
    if fit_NIM_flag
        % NIM (just playing here)
        nim0 = NIM( stim_par, repmat({'rectlin'},2,1), [1 -1], 'd2xt', 100 );
        k0 = nim0.subunits(1).filtK;
        nim0.subunits(1).filtK = k0;
        nim0.subunits(2).filtK = -k0;
        ksh = NIM.shift_mat_zpad(reshape(k0, [num_lags, NX*NY]), 3);
        nim0.subunits(2).filtK = reshape(ksh, [num_lags*NX*NY, 1]);
        nim0 = nim0.fit_filters( Rsh(:,rgc), Xstim, Ui, 'silent', silent );
        % nim0 = nim0.reg_path( Rsh(:,cc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' ); % to optimize regularization

        nim0b = nim0.fit_filters( Rsh(:,rgc), Xstim, Ui, 'fit_offsets', 1 );
        % Fitting offsets had good effect on main excitatory kernel

        % Add spike-history term
        nim1 = nim0b.init_spkhist( 12, 'doubling_time', 4, 'negcon');
        nim1 = nim1.fit_filters( Rsh(:,rgc), Xstim, Ui, 'fit_offsets', 1, 'silent', silent );
        nim1.display_model('Xstims', Xstim)
        drawnow

        fit_struct(rgc).nim_fit = nim1;
    end
  
end
toc

cd ~/Desktop/
save fit_struct fit_struct

%% Adding time-lagged spike trains from other cells (causal)

% loop over RGCs
for rgc = 1:NC
    disp(['%%%% FITTING RGC: ', num2str(rgc), ' w/ COUPLING TERMS %%%'])

    
    % just demonstrating GLM here
    other_ccs = setdiff(1:NC, rgc); % other neurons
    
    num_coupling_lags = 24;  % will be at lower resolution
    coupling_par = NIM.create_stim_params( [num_coupling_lags NC-1 1], 'stim_dt', refresh,'tent_spacing', frac);
    Cstim = NIM.create_time_embedding( Robs(:,other_ccs), coupling_par );
    % makes a time-lagged version of the other spike trains
    
    % Assemble design matrices
    Xs{1} = Xstim;
    Xs{2} =  Cstim;
    
    % make a coupled GLM model
    Pglm1 = fit_struct(rgc).glm_fit;
    Pglm1.stim_params(2) = coupling_par;  % add second stim params
    Pglm1 = Pglm1.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
    Pglm1 = Pglm1.fit_filters( Rsh(:,rgc), Xs, Ui );
    Pglm1.display_model('Xstims', Xs)
    drawnow
    fit_struct(rgc).coupled_glm = Pglm1;
    
    if fit_NIM_flag
        % make a couple NIM model (with auto-spike feedback)
        Pnim = fit_struct(rgc).nim_fit;
        Pnim.stim_params(2) = coupling_par;  % add second stim params
        Pnim = Pnim.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
        Pnim = Pnim.fit_filters( Rsh(:,rgc), Xs, Ui );
        Pnim.display_model('Xstims', Xs)
        fit_struct(rgc).coupled_nim = Pnim; 
        drawnow
    end
        
    % The big thing here is that a second stim-params was added, which corresponds to a different 'Xtarget'
    % Each subunit can be specified to act on a different stim based on its setting 'xtarg' (i.e. Pglm1.subunits(1).xtarg)
    % the model can be designed and fit from scratch as well, i.e.
    % >> Pglm1 = NIM( [stim_par coupling_par], {'lin', 'lin'}, [1 1], 'Xtargets', [1 2], ...)
   
end


cd ~/Desktop/
save fit_struct fit_struct


%% Model performance
[LLs,~,~,LLinfo] = glm0.eval_model(Rsh(:,1), Xstim, Xi ); LLnull = LLinfo.nullLL;
LLs(2) = glm1.eval_model(Rsh(:,1), Xstim, Xi );
LLs(3) = nim1.eval_model(Rsh(:,1), Xstim, Xi );
LLs(4) = Pglm1.eval_model(Rsh(:,1), Xs, Xi );
LLs(5) = Pnim.eval_model(Rsh(:,1), Xs, Xi );
LLs-LLnull


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECODING
[repX, repY, rep_frames] = size(repeat_movie);
repeat_movie(repeat_movie < 0.5) = -1;
repeat_movie(repeat_movie > 0.5) = 1;
rep_stim = reshape(repeat_movie, [repX * repY, rep_frames])';

% extract 5 seconds (1 trial) from repeats
tmp_rgc = 1; % choose an RGC
start_index = 20; % choose a repeat trial
% get spike_times
for rgc = 1:NC
    temp_spikes = repeat_spikes{rgc};
    decode_spikes = temp_spikes(temp_spikes >= starts(start_index) & temp_spikes < starts(start_index)+(rep_frames * refresh));
    dum = histc(decode_spikes, starts(start_index)+(0:(NTR-1))*(refresh/frac));
	decode_Robs(:,rgc) = dum;
end


rep_stim_par = NIM.create_stim_params( [num_lags repX repY], 'stim_dt', refresh, 'upsampling', frac, 'tent_spacing', frac);
decode_stim = NIM.create_time_embedding( rep_stim, rep_stim_par ); % design stim matrix
NTR = size(rep_stim,1) * frac;
%RobsR = NIM.shift_mat_zpad(	histc(decode_spikes, starts(start_index)+(0:(NTR-1))*(refresh/frac)), 0); % 

% Try evaluating the glm model

% get LN fit out of fit_struct
LN_fit = fit_struct(tmp_rgc).LN_fit;
%evaluate the LN fit given response
[LL, ~, ~, ~] = LN_fit.eval_model( decode_Robs(:,tmp_rgc), decode_stim);
LL

% get LN fit out of fit_struct
glm_fit = fit_struct(tmp_rgc).glm_fit;
%evaluate the LN fit given response
[LL, ~, ~, ~] = glm_fit.eval_model( decode_Robs(:,tmp_rgc), decode_stim);
LL

%% pop decoding

other_ccs = setdiff(1:NC, tmp_rgc); % other neurons
num_coupling_lags = 24;  % will be at lower resolution
coupling_par = NIM.create_stim_params( [num_coupling_lags NC-1 1], 'stim_dt', refresh,'tent_spacing', frac);
Cstim = NIM.create_time_embedding( decode_Robs(:,other_ccs), coupling_par );

% Assemble design matrices
Xs{1} = decode_stim;
Xs{2} =  Cstim;

% get LN fit out of fit_struct
coupled_glm = fit_struct(tmp_rgc).coupled_glm;
%evaluate the LN fit given response
[LL, ~, ~, ~] = coupled_glm.eval_model( RobsR, Xs);
LL









