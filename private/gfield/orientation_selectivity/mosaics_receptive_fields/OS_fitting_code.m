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

%% Inspect rasters

% define cell of interest
cell_id = 287;
[mapped_list, ~] = map_ei(grating_datarun, repeat_datarun, 'master_cell_type', cell_id, 'corr_threshold', .8);

temp_cell_index = get_cell_indices(grating_datarun, cell_id);
if isempty(mapped_list{temp_cell_index}); warning('cell did not map to white noise run'); end
wn_id = mapped_list{temp_cell_index};

wn_index = get_cell_indices(repeat_datarun, wn_id);
temp_spikes = repeat_datarun.spikes{wn_index};

reshaped_rep_trigs = reshape(rep_trigs, [], 1);
start_indices = [1:3:length(reshaped_rep_trigs)];

figure(2); clf;
temp_raster = get_raster(temp_spikes, reshaped_rep_trigs(start_indices), 'stop', 5, 'foa', 2);
axis tight

figure(3); clf;
temp_bin_size = rep_length / (300 * up_samp_fac);
[temp_psth, psth_bins] = get_psth(temp_spikes, reshaped_rep_trigs(start_indices), 'stop', 5, 'bin_size', temp_bin_size, 'foa', 3);
plot(psth_bins, temp_psth ./ max(temp_psth), 'b');

%% Compute the STA for a cell of interest

temp_cell_index = get_cell_indices(grating_datarun, cell_id);
if isempty(mapped_list{temp_cell_index}); warning('cell did not map to white noise run'); end
wn_id = mapped_list{temp_cell_index};
wn_index = get_cell_indices(repeat_datarun, wn_id);
temp_spikes = repeat_datarun.spikes{wn_index};


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


temp_sta = zeros(ni*nj, nLags);
for fm = nLags:1:size(Stim,2)
    if Robs(fm) > 0
        tmp_mov = Stim(:,fm-nLags+1:fm) .* Robs(fm);
        temp_sta = temp_sta + tmp_mov;
    end
end
STA = temp_sta ./ sum(Robs);

vis_sta = STA ./ abs(ext(STA(:)));
vis_sta = vis_sta + 1;
vis_sta = vis_sta ./ 2;

vis_sta = reshape(vis_sta, [ni, nj, nLags]);
vis_sta = repmat(vis_sta, [1 1 1 3]);
vis_sta = permute(vis_sta, [1,2,4,3]);
figure(1); clf;
for fm = 1:nLags
    image(squeeze(vis_sta(:,:,:,fm)))
    colormap gray
    axis square
    pause(0.1)
end

rSTA = reshape(STA, [ni, nj, nLags]);

%% get TC
stix_sta = repmat(rSTA, [1 1 1 3]);
% expand color dimension to make work with STA functions
stix_sta = permute(stix_sta, [1 2 4 3]);
sig_stix = significant_stixels(stix_sta, 'thresh', 3.5);
temp_tc = time_course_from_sta(stix_sta, sig_stix);
% colapse color dimension
temp_tc = temp_tc(:,1);
figure(1); clf;
plot(temp_tc)
axis square

%% Fit LNP and GLM Fitting

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


% Fit stimulus filter
disp('Fitting GLM: 1 excitatory filter')
L_fit = L_fit.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );
% look at filter
%k0 = reshape(L_fit.subunits(1).filtK, [nLags, NX*NY]);
%[mval,bestlag] = max(max(abs(k0)'));
L_fit.display_model('Xstims', Xstim)

% spike static nonlinearity
LN_fit = L_fit.fit_spkNL( Robs, Xstim, Ui );
LN_fit.display_model('Xstims', Xstim)

% Try spike history term
GLM_fit = LN_fit.init_spkhist(8, 'doubling_time',2, 'negcon');
GLM_fit = GLM_fit.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
GLM_fit.display_model('Xstims', Xstim)
GLM_fit = GLM_fit.fit_spkNL( Robs, Xstim, Ui );
GLM_fit.display_model('Xstims', Xstim)

cd ~/Desktop
exportgraphics(gcf, 'glm_model.pdf', 'ContentType', 'vector')

%% fit a more complex NIM

% initialize stuff
mod_signs = [1 1 1 -1]; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
NL_types = repmat({'rectlin'},1,length(mod_signs)); %make all upstream NLs threshold-linear
fit1 = NIM( params_stim, NL_types, mod_signs, 'd2xt', 100 ); 
optim_params.maxIter = 400;

% Fit stimulus filters
fit1 = fit1.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );
fit1.display_model('Xstims', Xstim)

% fit the spike (static) nonlinearity
fit2 = fit1.fit_spkNL( Robs, Xstim, Ui );
fit2.display_model('Xstims', Xstim)

% fit spike history
fit3 = fit2.init_spkhist(8, 'doubling_time',2, 'negcon');
fit3 = fit3.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
fit3.display_model('Xstims', Xstim)
% refit spike NL after fitting spike history
fit3 = fit3.fit_spkNL( Robs, Xstim, Ui );
fit3.display_model('Xstims', Xstim)

% fit subunit NLs
fit4 = fit3.init_nonpar_NLs( Xstim );
fit4 = fit4.fit_upstreamNLs( Robs, Xstim, Ui );
fit4.display_model('Xstims', Xstim)

exportgraphics(gcf, 'nim_model.pdf', 'ContentType', 'vector')

% refit spike history and static NL after fitting subunit NLs
%NIM_fit = fit4.init_spkhist(6, 'negcon');
NIM_fit = fit4.init_spkhist(8, 'doubling_time',2, 'negcon');
NIM_fit = NIM_fit.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
% refit spike NL after fitting spike history
NIM_fit = NIM_fit.fit_spkNL( Robs, Xstim, Ui );
NIM_fit.display_model('Xstims', Xstim)



%% PSTH Comparison: look at predicted firing rate from fit

model_to_examine1 = LN_fit;
model_to_examine2 = NIM_fit;

Xstim_rep = Xstim(1:length(psth_bins), :);
[LL, pred_rate, mod_internals, LL_data] = model_to_examine1.eval_model( [], Xstim_rep);
figure(3); clf; hold on
plot(psth_bins, temp_psth ./ max(temp_psth), 'k');
plot(psth_bins', pred_rate ./ max(pred_rate), 'r');
[LL, pred_rate, mod_internals, LL_data] = model_to_examine2.eval_model( [], Xstim_rep);
plot(psth_bins', pred_rate ./ max(pred_rate), 'b');
hold off

% compute R_squared
% test a circular shift to see if timing is off
c_shift = -1 * up_samp_fac;
norm_prediction = circshift(pred_rate, c_shift) ./ max(pred_rate);
norm_psth = temp_psth ./ max(temp_psth);
% cut off first nLags
sse = sum((norm_psth(nLags:end) - norm_prediction(nLags:end)').^2);
sst = sumsqr(norm_psth(nLags:end)) - mean(norm_psth(nLags:end)).^2;
r_squared = 1-(sse/sst)
 

corrcoef(norm_psth, norm_prediction)
temp_cc = corrcoef(norm_psth(nLags:end-nLags), norm_prediction(nLags:end-nLags))
gtext([num2str(temp_cc(1,2)),' , ',num2str(r_squared)]);
print(3, '~/Desktop/GLM_pred', '-dpdf')

rep_frame_num = size(Xstim_rep,1) * up_samp_fac;
rep_hist = zeros(length(start_indices), rep_frame_num);
for rp = 1:length(start_indices)
    
    bg_rep_time = reshaped_rep_trigs(start_indices(rp));
    rep_time_end = bg_rep_time + rep_length;
    bin_size = (rep_time_end - bg_rep_time) ./ rep_frame_num;

    temp_bins = bg_rep_time : bin_size : rep_time_end;
    temp_hist =histcounts(temp_spikes, temp_bins); 
    
    rep_hist(rp,:) = temp_hist;
end
%check_psth = sum(rep_hist);
%plot(psth_bins, check_psth ./ max(check_psth), 'g')

[LLs,LLnulls,pred_rates] = model_to_examine.eval_model_reps(rep_hist', Xstim_rep);



%% Raster Comparison
% generate raster from fit
num_sim_reps = 200;
[fit0nl_spiksR, fit0_psth] = nim_simulate(model_to_examine, Xstim_rep, num_sim_reps);
nim_trig_indices = find(fit0nl_spiksR == -1);
figure(2); subplot(2, 1, 2)
%RASTER
line_space = 0.2;
for raster_rps = 1:num_sim_reps
    if raster_rps == 1
        ith_trial_times = fit0nl_spiksR(1:nim_trig_indices(raster_rps)-1);
        line([ith_trial_times, ith_trial_times]', repmat([raster_rps-1+line_space, raster_rps-line_space], size(ith_trial_times,1), 1)',...
        'color', 'k', 'LineWidth', 1.0)
        %plot(fit0nl_spiksR(1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    else
        ith_trial_times = fit0nl_spiksR(nim_trig_indices(raster_rps-1)+1:nim_trig_indices(raster_rps)-1);
        line([ith_trial_times, ith_trial_times]', repmat([raster_rps-1+line_space, raster_rps-line_space], size(ith_trial_times,1), 1)',...
        'color', 'k', 'LineWidth', 1.0)
        %plot(fit0nl_spiksR(nim_trig_indices(raster_rps-1)+1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    end
    hold on
end
title('LNP')
%axis([0.5 5 0 num_sim_reps])
hold off
