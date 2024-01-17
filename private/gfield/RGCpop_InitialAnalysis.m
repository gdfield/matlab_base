%% Field Lab multi-neuron model fits -- initial demonstration
cd ~/Desktop/
load white_noise_example_data.mat

% Establish experimental parameters
NC = length(white_noise_spikes)
[NX, NY, NT] = size(white_noise_movie)
% organize stim to be +/- 1
stim = reshape(white_noise_movie,[NX*NY,NT])';
stim(stim < 0.5) = -1;
stim(stim > 0.5) = 1;

% Bin spike train
tbins = (0:NT)*refresh+white_noise_movie_start;
Robs = zeros(NT, NC);
% Binned at double resolution
tbins2 = (0:2*NT)*refresh/2+white_noise_movie_start;
Robs2 = zeros(2*NT, NC); % second version at double-resolution (8 ms)
% For 8x resolution
frac = 8;
Robs8 = zeros(frac*NT, NC); % second version at double-resolution (2 ms)
tbins8 = (0:frac*NT)*refresh/frac+white_noise_movie_start;

for nn = 1:10
	dum = histc(white_noise_spikes{nn}, tbins);
	Robs(:,nn) = dum(1:end-1);
	%dum = histc(white_noise_spikes{nn}, tbins2);
	%Robs2(:,nn) = dum(1:end-1);
	%dum = histc(white_noise_spikes{nn}, tbins8);
	%Robs8(:,nn) = dum(1:end-1);
end

% Save processed data
%save data_proc.mat stim NX NY Robs Robs2 Robs8 refresh

%% Model setup
lat_skip = 2;  % how many time lags to throw out preceding stimulus
num_lags = 24; % number of time lags in temporal kernel: note that response is analyzed here at refresh/2

% for no upsampling
stim_par = NIM.create_stim_params( [num_lags NX NY], 'stim_dt', refresh, 'upsampling', 1 );
Xstim = NIM.create_time_embedding( stim, stim_par ); % design stim matrix
Rsh = NIM.shift_mat_zpad( Robs, -lat_skip );  % shifted responses
[Ui, Xi] = NIM.generate_XVfolds( NT ); % establish 5-fold cross-validation indices

% for upsampling 2x
stim_par = NIM.create_stim_params( [num_lags NX NY], 'stim_dt', refresh, 'upsampling', 2 );
Xstim = NIM.create_time_embedding( stim, stim_par ); % design stim matrix
Rsh = NIM.shift_mat_zpad( Robs2, -lat_skip );  % shifted responses
[Ui, Xi] = NIM.generate_XVfolds( NT*2 ); % establish 5-fold cross-validation indices



% STAs
stas = Xstim'*Rsh;
figure
for cc = 1:10
	subplot(3,4, cc)
	k = reshape(stas(:,cc), [num_lags, NX*NY]);
	plot(k)
end

% GLM
cc = 1;  % fit example cell

glm0 = NIM( stim_par, {'lin'}, 1, 'd2xt', 100 );
glm0 = glm0.fit_filters( Rsh(:,cc), Xstim, Ui );
glm0 = glm0.reg_path( Rsh(:,cc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' );

% Add spike-history term
glm1 = glm0.init_spkhist( 12, 'doubling_time', 4);
glm1 = glm1.fit_filters( Rsh(:,cc), Xstim, Ui );

% NIM (just playing here)
nim0 = NIM( stim_par, repmat({'rectlin'},2,1), [1 -1], 'd2xt', 100 );
k0 = nim0.subunits(1).filtK;
nim0.subunits(1).filtK = k0;
nim0.subunits(2).filtK = -k0;
ksh = NIM.shift_mat_zpad(reshape(k0, [num_lags, NX*NY]), 3);
nim0.subunits(2).filtK = reshape(ksh, [num_lags*NX*NY, 1]);
nim0 = nim0.fit_filters( Rsh(:,cc), Xstim, Ui );
% nim0 = nim0.reg_path( Rsh(:,cc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' ); % to optimize regularization

nim0b = nim0.fit_filters( Rsh(:,cc), Xstim, Ui, 'fit_offsets', 1 );
nim0b.display_model('Xstims', Xstim)
% Fitting offsets had good effect on main excitatory kernel

% Add spike-history term
nim1 = nim0b.init_spkhist( 12, 'doubling_time', 4, 'negcon');
nim1 = nim1.fit_filters( Rsh(:,cc), Xstim, Ui, 'fit_offsets', 1 );

%% Adding time-lagged spike trains from other cells (causal)
% just demonstrating GLM here
other_ccs = setdiff(1:NC, cc); % other neurons

num_coupling_lags = 24;  % will be at lower resolution
coupling_par = NIM.create_stim_params( [num_coupling_lags NC-1 1], 'stim_dt', refresh, 'tent_spacing', 2 );
Cstim = NIM.create_time_embedding( Robs(:,other_ccs), coupling_par );
% makes a time-lagged version of the other spike trains

% Assemble design matrices
Xs{1} = Xstim;
Xs{2} =  Cstim;

% make a coupled GLM model
Pglm1 = glm1; 
Pglm1.stim_params(2) = coupling_par;  % add second stim params
Pglm1 = Pglm1.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
Pglm1 = Pglm1.fit_filters( Rsh(:,cc), Xs, Ui );
Pglm1.display_model('Xstims', Xs)

% make a couple NIM model (with auto-spike feedback)
Pnim = nim1; 
Pnim.stim_params(2) = coupling_par;  % add second stim params
Pnim = Pnim.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
Pnim = Pnim.fit_filters( Rsh(:,cc), Xs, Ui );
Pnim.display_model('Xstims', Xs)

% The big thing here is that a second stim-params was added, which corresponds to a different 'Xtarget'
% Each subunit can be specified to act on a different stim based on its setting 'xtarg' (i.e. Pglm1.subunits(1).xtarg)
% the model can be designed and fit from scratch as well, i.e.
% >> Pglm1 = NIM( [stim_par coupling_par], {'lin', 'lin'}, [1 1], 'Xtargets', [1 2], ...)

%% Model performance
[LLs,~,~,LLinfo] = glm0.eval_model(Rsh(:,cc), Xstim, Xi ); LLnull = LLinfo.nullLL;
LLs(2) = glm1.eval_model(Rsh(:,cc), Xstim, Xi );
LLs(3) = nim1.eval_model(Rsh(:,cc), Xstim, Xi );
LLs(4) = Pglm1.eval_model(Rsh(:,cc), Xs, Xi );
LLs(5) = Pnim.eval_model(Rsh(:,cc), Xs, Xi );
LLs-LLnull
