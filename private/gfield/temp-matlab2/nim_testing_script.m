%% Load data

%datapath = '/Volumes/Disk2/Data/2019-06-24-0/data005_MR/data005_MR';
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-06-24-0/data005_MR/data005_MR';
moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
%moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

%% Choose neuron, and plot raster
cell_id = 4877; % example OFF brisk sustained
% cell_id = 1606; % OFF brisk sustained
% cell_id = 438; % ON brisk sustained
% cell_id = 781; % example ON brisk transient

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, cell_id, 'marks_params', marks_params);
cell_index = get_cell_indices(datarun, cell_id);
temp_spikes = datarun.spikes{cell_index};

% get center of mass and of STA and cut out an ROI to reduct dim of stim
temp_com = abs([0 40] - round(datarun.vision.sta_fits{cell_index}.mean)) +[1 1];
temp_size = 3*sqrt(prod(datarun.vision.sta_fits{cell_index}.sd));
window_size = round(temp_size);
bg_height = floor(temp_com(2) - window_size);
end_height = floor(temp_com(2) + window_size);
bg_width = floor(temp_com(1) - window_size);
end_width = floor(temp_com(1) + window_size);
wind_dim = window_size*2 + 1;

frame_num = 70000;
[mov,height,width,duration,refresh] = get_movie(moviepath, datarun.triggers, frame_num);
mov = squeeze(mov(bg_height:end_height, bg_width:end_width,1,:));
mov = reshape(mov, [(wind_dim * wind_dim), frame_num]);
mov = permute(mov, [2 1]) -0.5;

%mov = squeeze(mov(:, :, 1, :));
%mov = reshape(mov, [height * width, frame_num]);
%mov = permute(mov, [2 1]) -0.5;
dt = refresh./ 1000;



frame_time = refresh * frame_num /1000;
cut_spikes = temp_spikes(temp_spikes < frame_time);

%% Compute the NL
%datarun = get_snls(datarun, cell_id, 'fit', 'exp', 'new', true);
%SNL_fit = datarun.stas.snls{cell_index}.fit_params;

%% Fit GLM model
NIM_params.nLags = 24; %number of time lags for estimating stimulus filters
NIM_params.dt = refresh./ 1000;
NIM_params.up_samp_fac = 1; % temporal up-sampling factor applied to stimulus 
NIM_params.tent_basis_spacing = 1; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
NIM_params.silent = 0;
NIM_params.binsize = NIM_params.dt/NIM_params.up_samp_fac;
NIM_params.model_to_fit = [1 -1];
%NIM_params.wind_dim = wind_dim;

%NIM_params.NX = 53;
%NIM_params.NY = 40;

NIM_params.NX = wind_dim;
NIM_params.NY = wind_dim;

fit_struct = NIM_wrapper_glm_only(mov, cut_spikes, NIM_params);

%% Load DATA for WN repeats

datapath = '/Users/gfield/Analysis/2017-01-16-0/data005-map_KR/data005-map_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun_rep = load_data(datapath);
datarun_rep = load_neurons(datarun_rep);
datarun_rep = load_params(datarun_rep);
datarun_rep = load_java_movie(datarun_rep, moviepath);

% get cell index for cell_id specified above: make sure mapping worked!!!
cell_index_rep = get_cell_indices(datarun_rep, cell_id);


%% Plot raster and PSTH for DATA
figure(4); clf
subplot(4,1,1)
temp_raster = get_raster(datarun_rep.spikes{cell_index_rep}, datarun_rep.triggers([7:3:597]), 'plot', false);
%RASTER
for raster_rps = 1:size(temp_raster, 1)
    plot(temp_raster{raster_rps},raster_rps, '.k')
    hold on
end
title('data')
axis([0.5 5 0 200])
hold off
%PSTH
[tmp_psth, psth_bins] = get_psth(datarun_rep.spikes{cell_index_rep}, datarun_rep.triggers([1:3:597]), 'bin_size', (1/60.35), 'stop', 5.0);
figure(5); clf; hold on
plot(tmp_psth(1:302), 'k')
axis tight


%% Extract stimulus sequence for WN repeats
repeat_length = 5; % units of seconds
monitor_refresh = 60.35; %Hz
% get center of mass and of STA and cut out an ROI to reduce dim of stim
frame_num = floor(monitor_refresh*repeat_length);
[mov_rpts,height,width,duration,refresh] = get_movie(moviepath, datarun_rep.triggers(1:3), frame_num);
mov_rpts = squeeze(mov_rpts(bg_height:end_height, bg_width:end_width,1,:));
mov_rpts = reshape(mov_rpts, [(wind_dim * wind_dim), frame_num]);
mov_rpts = permute(mov_rpts, [2 1]) -0.5;


params_stim = NIM.create_stim_params([NIM_params.nLags NIM_params.wind_dim, NIM_params.wind_dim],...
    'stim_dt',	NIM_params.dt, 'upsampling', NIM_params.up_samp_fac, 'tent_spacing', NIM_params.tent_basis_spacing);


Nreps = length(temp_raster);
% make raster for NIM
spks_obs = [];
for rr = 1:Nreps
	spks_obs = [spks_obs; temp_raster{rr}; -1;];
end
NTR = size(mov_rpts,1) * NIM_params.up_samp_fac;
% Bin to make Robs -- note shift by latenct change
RobsR = NIM.shift_mat_zpad(	histc(spks_obs,(0:(NTR-1))*NIM_params.binsize)/Nreps, 2*NIM_params.up_samp_fac);
Xstim_rep = NIM.create_time_embedding(mov_rpts,params_stim); 


%%

RobsR_norm = RobsR ./ NIM_params.dt;
figure(5); clf;
plot(RobsR_norm, 'k')
hold on

dwnsample_RobsR = zeros(length(RobsR)/5,1);
for cc = 1:length(RobsR)/5;
    bg = (5*(cc-1)) + 1;
    nd = (5*cc);
    dwnsample_RobsR(cc) = sum(RobsR(bg:nd));
end
figure(6); clf;
plot(dwnsample_RobsR./ NIM_params.dt, 'k')
hold on

%% create raster & PSTH for fit0;
num_sim_reps = 100;
[fit0_spiksR, fit0_psth] = nim_simulate(fit_struct.fit0, Xstim_rep, num_sim_reps);
nim_trig_indices = find(fit0_spiksR == -1);
figure(4)
subplot(4,1,2)
%RASTER
for raster_rps = 1:num_sim_reps
    if raster_rps == 1
        plot(fit0_spiksR(1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    else
        plot(fit0_spiksR(nim_trig_indices(raster_rps-1)+1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    end
    hold on
end
title('LNP')
axis([0.5 5 0 num_sim_reps])
hold off

figure(5);
plot(fit0_psth./ NIM_params.dt, 'b')

dwnsample_fit0_psth = zeros(length(RobsR)/5,1);
for cc = 1:length(RobsR)/5
    bg = (5*(cc-1)) + 1;
    nd = (5*cc);
    dwnsample_fit0_psth(cc) = sum(fit0_psth(bg:nd));
end
figure(6)
plot(dwnsample_fit0_psth./ NIM_params.dt, 'b')


%% ccreate raster & PSTH for fitS;
num_sim_reps = 100;
[fitS_spiksR, fitS_psth] = nim_simulate(fit_struct.fitS, Xstim_rep, num_sim_reps);
nim_trig_indices = find(fitS_spiksR == -1);
figure(4)
subplot(4,1,3)
%RASTER
for raster_rps = 1:num_sim_reps
    if raster_rps == 1
        plot(fitS_spiksR(1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    else
        plot(fitS_spiksR(nim_trig_indices(raster_rps-1)+1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    end
    hold on
end
title('GLM')
axis([0.5 5 0 num_sim_reps])
hold off

figure(5);
plot(fitS_psth./ NIM_params.dt, 'g')

dwnsample_fitS_psth = zeros(length(RobsR)/5,1);
for cc = 1:length(RobsR)/5
    bg = (5*(cc-1)) + 1;
    nd = (5*cc);
    dwnsample_fitS_psth(cc) = sum(fitS_psth(bg:nd));
end
figure(6)
plot(dwnsample_fitS_psth./ NIM_params.dt, 'g')


%% create raster & PSTH for fit1S;
num_sim_reps = 100;
[fit1S_spiksR, fit1S_psth] = nim_simulate(fit_struct.fit1S, Xstim_rep, num_sim_reps);
nim_trig_indices = find(fit1S_spiksR == -1);
figure(4)
subplot(4,1,4)
%RASTER
for raster_rps = 1:num_sim_reps
    if raster_rps == 1
        plot(fit1S_spiksR(1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    else
        plot(fit1S_spiksR(nim_trig_indices(raster_rps-1)+1:nim_trig_indices(raster_rps)-1),raster_rps, '.k')
    end
    hold on
end
title('NIM [1 -1]')
axis([0.5 5 0 num_sim_reps])
hold off

figure(5);
plot(fit1S_psth./ NIM_params.dt, 'r')

dwnsample_fit1S_psth = zeros(length(RobsR)/5,1);
for cc = 1:length(RobsR)/5
    bg = (5*(cc-1)) + 1;
    nd = (5*cc);
    dwnsample_fit1S_psth(cc) = sum(fit1S_psth(bg:nd));
end
figure(6)
plot(dwnsample_fit1S_psth./ NIM_params.dt, 'r')
legend('data', 'LNP', 'GLM', 'NIM 1,-1')


%%

% Dan's code for comparison:
R2s = 1-[var(RobsR-fit0_psth) var(RobsR-fit1S_psth)]/var(RobsR)

% Kiersten's code for comparison
perf = get_pred_fr_metrics(RobsR', [fit0_psth, fitS_psth fit1S_psth]', [1 1 1 1])

perf = get_pred_fr_metrics(dwnsample_RobsR', [dwnsample_fit0_psth, dwnsample_fitS_psth dwnsample_fit1S_psth]', [1 1 1 1])





