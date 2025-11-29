%% load data
grating_path = '/Users/gfield/Analysis/2018-11-30-0/data004_KR-map/data004_KR-map';
grating_datarun = load_data(grating_path);
grating_datarun = load_neurons(grating_datarun);
grating_datarun = load_params(grating_datarun);
grating_datarun = load_ei(grating_datarun, 'all');

%load list of os cells in this run
cd ~/Desktop/2018-11-30/
load os_cell_list

% load WN repeat run
repeat_path = '/Users/gfield/Analysis/2018-11-30-0/data002_KR-map_KR/data002_KR-map_KR';
repeat_datarun = load_data(repeat_path);
repeat_datarun = load_neurons(repeat_datarun);
repeat_datarun = load_params(repeat_datarun);
repeat_datarun = load_ei(repeat_datarun, 'all');

% map from grainting run to the white noise + repeats run
[mapped_list, ~] = map_ei(grating_datarun, repeat_datarun, 'master_cell_type', os_cell_list, 'corr_threshold', .9);

% organize interleaves
%from stimulus file:
total_rep_time = 1125; %in seconds
total_nonrep_time = 4050;
%repeat lengths set to 125s
num_interleaves = total_rep_time/125; %should evenly divide
nonrep_length = total_nonrep_time/num_interleaves;
rep_length = 5; %in seconds, one repeat
num_reps_win = 125/rep_length; %evenly divide

% some fitting parameters
nLags = 24;
up_samp_fac = 1;
tent_basis_spacing = 1;
dt = 1000/60.35;

%% Sort out the triggers and the movies
triggers = repeat_datarun.triggers;
non_rep_trigs = zeros(270, num_interleaves); %450*60/100 per interleave     1:(75+270):2070
rep_trigs = zeros(75, num_interleaves); %125*60/100 per interleave

for ni = 1:num_interleaves
    non_rep_trigs(:,ni) = triggers( 345*(ni-1)+1 : 345*(ni-1)+1+269 );
    rep_trigs(:,ni) = triggers( 345*(ni-1)+1+270 : 345*(ni-1)+1+270+74 );
end     

seeds = [11111 22222 33333 44444 55555 66666 77777 88888 99999];
nm1s = '/Users/gfield/Development/movie-xml2/BW-60-1-0.48-';
nm3s = '-13x10-60.35_xoffset.xml';

tr = str2num(nm1s(end-6));
for ni = 1:num_interleaves
    nm2 = seeds(ni);
    [temp_mov,~,~,dur,frame_refresh] = get_movie([nm1s num2str(nm2) nm3s],non_rep_trigs(:,ni),nonrep_length*60/tr); %num of unique frames
    mov{ni} = squeeze(temp_mov(:,:,1,:));
end

ni = size(temp_mov,1);
nj = size(temp_mov,2);
white_noise_movie = cat(3,mov{1},mov{2},mov{3},mov{4},mov{5},mov{6},mov{7},mov{8},mov{9}); %do this smarter
Stim0 = reshape(white_noise_movie,ni*nj,[]);
Stim = bsxfun(@minus,Stim0,0.5);

%%
cell_list = [3001];
num_rgcs = length(cell_list);

model_types = {[1 1], [1 -1], [1 1 1], [1 1 -1], [1 -1 -1], [1 1 1 1], [1 1 1 -1], [1 1 -1 -1]};

fit_structures = struct([]);

for rgc = 1:num_rgcs
    
    disp(['******** Fitting cell ', num2str(cell_list(rgc)), ' *********'])

    %% Compute the STA for a cell of interest

    cell_id = cell_list(rgc);
    [mapped_list, ~] = map_ei(grating_datarun, repeat_datarun, 'master_cell_type', cell_id, 'corr_threshold', .9);
    temp_cell_index = get_cell_indices(grating_datarun, cell_id);
    if isempty(mapped_list{temp_cell_index}); warning('cell did not map to white noise run'); end
    wn_id = mapped_list{temp_cell_index};
    wn_index = get_cell_indices(repeat_datarun, wn_id);
    temp_spikes = repeat_datarun.spikes{wn_index};
    reshaped_rep_trigs = reshape(rep_trigs, [], 1);
    start_indices = [1:3:length(reshaped_rep_trigs)];

    temp_bin_size = rep_length / (300 * up_samp_fac);
    [temp_psth, psth_bins] = get_psth(temp_spikes, reshaped_rep_trigs(start_indices), 'stop', 5, 'bin_size', temp_bin_size, 'foa', 3);
    norm_psth = temp_psth ./ max(temp_psth);   
    %figure(3); clf;
    %plot(psth_bins, norm_psth, 'b');

    Robs = [];
    for n = 1:num_interleaves
        % expand the triggers to approximate frame times
        frame_bins = [];
        diff_trig_times = diff(non_rep_trigs(:,n));
        for j = 1:length(diff_trig_times)
            temp_vec = non_rep_trigs(j,n) : diff_trig_times(j)/(100*up_samp_fac  - 1) : non_rep_trigs(j+1,n);
            frame_bins = [frame_bins temp_vec];
        end
        %tag on another 100 frames for last trigger
        temp_vec = non_rep_trigs(j+1,n) : mean(diff_trig_times)/(100*up_samp_fac  - 1) : non_rep_trigs(j+1,n) + mean(diff_trig_times) + mean(diff_trig_times)/(100*up_samp_fac  - 1);
        frame_bins = [frame_bins temp_vec];

        [temp_hist, temp_bins] = histcounts(temp_spikes, frame_bins);

        Robs = [Robs, temp_hist];    
    end

    %% Fit LNP and GLM 

    % initialize some stuff...
    NX = ni;
    NY = nj;
    % Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
    [NFRAMES, nXPix] = size(Stim');
    NT = NFRAMES*up_samp_fac;
    % Init stim parameters structure
    params_stim = NIM.create_stim_params([nLags NX NY], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);
    % generate stimulus matrix
    Xstim = NIM.create_time_embedding(Stim',params_stim); 
    [Ui, XVi] = NIM.generate_XVfolds( NT );
    mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; % define subunit as linear 
    % New object-based definition of model (note regularization can be set here)
    L_fit = NIM( params_stim, NL_types, mod_signs, 'd2xt', 100 ); 
    % also using spatial regularization too -- can do seprately or as one
    optim_params.maxIter = 500;
    Xstim_rep = Xstim(1:length(psth_bins), :);

    % Fit stimulus filter
    disp('Fitting GLM: 1 excitatory filter')
    L_fit = L_fit.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );
    [LL, pred_rate, mod_internals, LL_data] = L_fit.eval_model( [], Xstim_rep);
    L_fit_struct.LL = LL; L_fit_struct.pred_rate = pred_rate; L_fL_fit_structit.LL_data = LL_data;
    norm_prediction = pred_rate ./ max(pred_rate);
    temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
    L_fit_struct.cor_coef = temp_corr(2);
   
    % spike static nonlinearity
    LN_fit = L_fit.fit_spkNL( Robs, Xstim, Ui );
    [LL, pred_rate, mod_internals, LL_data] = LN_fit.eval_model( [], Xstim_rep);
    LN_fit_struct.LL = LL; LN_fit_struct.pred_rate = pred_rate; LN_fit_struct.LL_data = LL_data;
    norm_prediction = pred_rate ./ max(pred_rate);
    temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
    LN_fit_struct.cor_coef = temp_corr(2);

    % Try spike history term
    GLM_fit = LN_fit.init_spkhist(8, 'doubling_time',2, 'negcon');
    GLM_fit = GLM_fit.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
    GLM_fit = GLM_fit.fit_spkNL( Robs, Xstim, Ui );
    [LL, pred_rate, mod_internals, LL_data] = GLM_fit.eval_model( [], Xstim_rep);
    GLM_fit_struct.LL = LL; GLM_fit_struct.pred_rate = pred_rate; GLM_fit_struct.LL_data = LL_data;
    norm_prediction = pred_rate ./ max(pred_rate);
    temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
    GLM_fit_struct.cor_coef = temp_corr(2);

    fit_structures(rgc).L_fit = L_fit;
    fit_structures(rgc).L_fit_struct = L_fit_struct;
    fit_structures(rgc).LN_fitt = LN_fit;
    fit_structures(rgc).LN_fit_struct = LN_fit_struct;
    fit_structures(rgc).GLM_fit = GLM_fit;
    fit_structures(rgc).GLM_fit_struct = GLM_fit_struct;

    % fit multi-filter NIMs
    for mt = 1:length(model_types)
        %% fit a more complex NIM
        % initialize stuff
        mod_signs = model_types{mt}; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
        NL_types = repmat({'rectlin'},1,length(mod_signs)); %make all upstream NLs threshold-linear
        fit1 = NIM( params_stim, NL_types, mod_signs, 'd2xt', 100 ); 
  
        % Fit stimulus filters
        fit1 = fit1.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );
        [LL, pred_rate, mod_internals, LL_data] = fit1.eval_model( [], Xstim_rep);
        fit1_struct.LL = LL; fit1_struct.pred_rate = pred_rate; fit1_struct.LL_data = LL_data;
        norm_prediction = pred_rate ./ max(pred_rate);
        temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
        fit1_struct.cor_coef = temp_corr(2);

        % fit the spike (static) nonlinearity
        fit2 = fit1.fit_spkNL( Robs, Xstim, Ui );
        [LL, pred_rate, mod_internals, LL_data] = fit2.eval_model( [], Xstim_rep);
        fit2_struct.LL = LL; fit2_struct.pred_rate = pred_rate; fit2_struct.LL_data = LL_data;
        norm_prediction = pred_rate ./ max(pred_rate);
        temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
        fit2_struct.cor_coef = temp_corr(2);

        % fit spike history
        fit3 = fit2.init_spkhist(8, 'doubling_time',2, 'negcon');
        fit3 = fit3.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
        % refit spike NL after fitting spike history
        fit3 = fit3.fit_spkNL( Robs, Xstim, Ui );
        [LL, pred_rate, mod_internals, LL_data] = fit3.eval_model( [], Xstim_rep);
        fit3_struct.LL = LL; fit3_struct.pred_rate = pred_rate; fit3_struct.LL_data = LL_data;
        norm_prediction = pred_rate ./ max(pred_rate);
        temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
        fit3_struct.cor_coef = temp_corr(2);

        % fit subunit NLs
        fit4 = fit3.init_nonpar_NLs( Xstim );
        fit4 = fit4.fit_upstreamNLs( Robs, Xstim, Ui );
        [LL, pred_rate, mod_internals, LL_data] = fit4.eval_model( [], Xstim_rep);
        fit4_struct.LL = LL; fit4_struct.pred_rate = pred_rate; fifit4_structt4.LL_data = LL_data;
        norm_prediction = pred_rate ./ max(pred_rate);
        temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
        fit4_struct.cor_coef = temp_corr(2);
        
        % refit spike history and static NL after fitting subunit NLs
        %NIM_fit = fit4.init_spkhist(6, 'negcon');
        NIM_fit = fit4.init_spkhist(8, 'doubling_time',2, 'negcon');
        NIM_fit = NIM_fit.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
        % refit spike NL after fitting spike history
        NIM_fit = NIM_fit.fit_spkNL( Robs, Xstim, Ui );
        [LL, pred_rate, mod_internals, LL_data] = NIM_fit.eval_model( [], Xstim_rep);
        NIM_fit_struct.LL = LL; NIM_fit_struct.pred_rate = pred_rate; NIM_fit_struct.LL_data = LL_data;
        norm_prediction = pred_rate ./ max(pred_rate);
        temp_corr = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags));
        NIM_fit_struct.cor_coef = temp_corr(2);

        fit_structures(rgc).NIM(mt).fit1 = fit1;
        fit_structures(rgc).NIM(mt).fit1_struct = fit1_struct;
        fit_structures(rgc).NIM(mt).fit2 = fit2;
        fit_structures(rgc).NIM(mt).fit2_struct = fit2_struct;
        fit_structures(rgc).NIM(mt).fit3 = fit3;
        fit_structures(rgc).NIM(mt).fit3_struct = fit3_struct;
        fit_structures(rgc).NIM(mt).fit4 = fit4;
        fit_structures(rgc).NIM(mt).fit4_struct = fit4_struct;
        fit_structures(rgc).NIM(mt).NIM_fit = NIM_fit;
        fit_structures(rgc).NIM(mt).NIM_fit_struct = NIM_fit_struct;

    end
end

