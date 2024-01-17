
%% load WN run
datapath = '/Users/gfield/Analysis/2017-01-16-0/data006_KR/data006_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_sta(datarun, 'load_sta', 'all');

%% load WN repeats run
datapath = '/Users/gfield/Analysis/2017-01-16-0/data005-map_KR/data005-map_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun_rep = load_data(datapath);
datarun_rep = load_neurons(datarun_rep);
datarun_rep = load_params(datarun_rep);
datarun_rep = load_ei(datarun_rep, 'all');
datarun_rep = load_java_movie(datarun_rep, moviepath);

repeat_length = 5; % units of seconds
monitor_refresh = 60.35; %Hz
frame_num = floor(monitor_refresh*repeat_length);
[mov_rpts,height,width,duration,refresh] = get_movie(moviepath, datarun_rep.triggers(1:3), frame_num);


%% load fits
load NIM_fits 

NIM_params.nLags = 25; %number of time lags for estimating stimulus filters
NIM_params.dt = refresh./ 1000;
NIM_params.up_samp_fac = 5; % temporal up-sampling factor applied to stimulus 
NIM_params.tent_basis_spacing = 5; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
NIM_params.silent = 0;
NIM_params.binsize = NIM_params.dt/NIM_params.up_samp_fac;
NIM_params.model_to_fit = [1 -1];


%% map cells
temp_type = {'off brisk transient'};
temp_indices = get_cell_indices(datarun, temp_type);
num_rgcs = length(temp_indices);

[cell_list, ~] = map_ei(datarun, datarun_rep, 'master_cell_type', temp_type);

% get best matches in WN repeats
conf_cell_list = [];
for rgc = 1:num_rgcs
    if ~isempty(cell_list{temp_indices(rgc)})
        conf_cell_list = [conf_cell_list, cell_list{temp_indices(rgc)}];
    else
        conf_cell_list = [conf_cell_list -1];
    end
end


model_perf = zeros(4,num_rgcs);
for rgc = 1:num_rgcs
    
    % ensure that the cell mapped, skip if not
    if conf_cell_list(rgc) == -1
        continue
    end
        
    %%% cut out window of movie %%%
    temp_com = abs([0 40] - round(datarun.vision.sta_fits{temp_indices(rgc)}.mean)) +[1 1];
    temp_size = 2.5*sqrt(prod(datarun.vision.sta_fits{temp_indices(rgc)}.sd));
    window_size = round(temp_size);
    bg_height = floor(temp_com(2) - window_size);
    end_height = floor(temp_com(2) + window_size);
    bg_width = floor(temp_com(1) - window_size);
    end_width = floor(temp_com(1) + window_size);    
    if bg_height < 1
        bg_height = 1;
    end
    if bg_width < 1
        bg_width = 1;
    end
    if end_height > height
        end_height = height;
    end
    if end_width > width
        end_width = width;
    end   
    NIM_params.NX = end_width - bg_width +1;
    NIM_params.NY = end_height - bg_height +1;
    mov_windowed = squeeze(mov_rpts(bg_height:end_height, bg_width:end_width,1,:));
    mov_windowed = reshape(mov_windowed, [(NIM_params.NX * NIM_params.NY), frame_num]);
    mov_windowed = permute(mov_windowed, [2 1]) -0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
    tmp_fit_struct = NIM_fits{temp_indices(rgc)};

    rep_index = get_cell_indices(datarun_rep, conf_cell_list(rgc));
    temp_raster = get_raster(datarun_rep.spikes{rep_index}, datarun_rep.triggers([7:3:597]), 'plot', false);
    
    %%% cut get the psth for the cell %%%
    params_stim = NIM.create_stim_params([NIM_params.nLags NIM_params.NX, NIM_params.NY],...
        'stim_dt',	NIM_params.dt, 'upsampling', NIM_params.up_samp_fac, 'tent_spacing', NIM_params.tent_basis_spacing);

    Nreps = length(temp_raster);
    % make raster for NIM
    spks_obs = [];
    for rr = 1:Nreps
        spks_obs = [spks_obs; temp_raster{rr}; -1;];
    end
    NTR = size(mov_windowed,1) * NIM_params.up_samp_fac;
    % Bin to make Robs -- note shift by latenct change
    RobsR = NIM.shift_mat_zpad(	histc(spks_obs,(0:(NTR-1))*NIM_params.binsize)/Nreps, 3*NIM_params.up_samp_fac);
    Xstim_rep = NIM.create_time_embedding(mov_windowed,params_stim); 

    %%% simulate different models %%%
    num_sim_reps = 200;
    [~, LNP_psth] = nim_simulate(tmp_fit_struct.fit0, Xstim_rep, num_sim_reps);
    [~, GLM_psth] = nim_simulate(tmp_fit_struct.fitS, Xstim_rep, num_sim_reps);
    [~, NIM_psth] = nim_simulate(tmp_fit_struct.fit1, Xstim_rep, num_sim_reps);
    [~, NIMH_psth] = nim_simulate(tmp_fit_struct.fit1S, Xstim_rep, num_sim_reps);

    if 1
        figure(1); clf;
        plot(RobsR, 'k')
        hold on
        plot(LNP_psth, 'r')
        plot(NIMH_psth, 'g')
        hold off
    end

    temp_perf = get_pred_fr_metrics(RobsR(150:1500)', [LNP_psth(150:1500), GLM_psth(150:1500) NIM_psth(150:1500) NIMH_psth(150:1500)]', [1 1 1 1]);   
 
    model_perf(1:4,rgc) = temp_perf(1:4,3);
end

figure(2); clf;
plot(model_perf, 'k.')
hold on
errorbar(median(model_perf'), std(model_perf')./sqrt(num_rgcs-1), 'sr')
axis([0 5 0 1])
title(temp_type{1})
ylabel('variance explained')
xlabel('model type: LNP, GLM, NIM, NIM+SH')


    
    
    
    
    
    
    