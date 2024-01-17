datapath = '/Users/gfield/Analysis/2017-01-16-0/data006_KR/data006_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_java_movie(datarun, moviepath);
datarun = load_ei(datarun, 'all');

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);

frame_num = 215000;
[mov,height,width,duration,refresh] = get_movie(moviepath, datarun.triggers, frame_num);
frame_time = refresh * frame_num /1000;

%%
load NIM_fits

%% OFF brisk transient cells
cell_type = {'ON brisk transient'};
temp_indices = get_cell_indices(datarun, cell_type);


nLags = 25; %number of time lags for estimating stimulus filters
dt = refresh./ 1000;
up_samp_fac = 2; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 1; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
silent = 0;
binsize = dt/up_samp_fac;



for rgc = 1:length(cell_ids);
    
    temp_NIM = NIM_fits{temp_indices(rgc)};
    
    % get center of mass and of STA and cut out an ROI to reduct dim of stim
    temp_com = abs([0 40] - round(datarun.vision.sta_fits{temp_indices(rgc)}.mean)) +[1 1];
    temp_size = 2.5*sqrt(prod(datarun.vision.sta_fits{temp_indices(rgc)}.sd));
    window_size = round(temp_size);
    bg_height = floor(temp_com(2) - window_size);
    end_height = floor(temp_com(2) + window_size);
    bg_width = floor(temp_com(1) - window_size);
    end_width = floor(temp_com(1) + window_size);
    wind_dim = window_size*2 + 1;

    mov_windowed = squeeze(mov(bg_height:end_height, bg_width:end_width,1,:));
    mov_windowed = reshape(mov_windowed, [(wind_dim * wind_dim), frame_num]);
    mov_windowed = permute(mov_windowed, [2 1]) -0.5;

    temp_spikes = datarun.spikes{temp_indices(rgc)};
    cut_spikes = temp_spikes(temp_spikes < frame_time);

    
    NX = wind_dim;

    % Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
    [NFRAMES, nXPix] = size(mov);
    NT = NFRAMES*up_samp_fac;

    % Init stim parameters structure
    params_stim = NIM.create_stim_params([nLags NX, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);

    % generate stimulus matrix
    Xstim = NIM.create_time_embedding(mov_windowed,params_stim); 

    % bin spike times 
    Robs = NIM.Spks2Robs(cut_spikes, binsize, NT );

    [Ui, XVi] = NIM.generate_XVfolds( NT );    
    
    
    
    % internal LLx (cross-validated likelihood, although was used for meta-params)
    [LLs,~,~,fp] = temp_NIM.fit0.eval_model(Robs, Xstim, XVi );
    LLs(2) = temp_NIM.fitS.eval_model(Robs, Xstim, XVi );
    LLs(3) = temp_NIM.fit1.eval_model(Robs, Xstim, XVi );
    LLs(4) = temp_NIM.fit1S.eval_model(Robs, Xstim, XVi );
    LLs-fp.nullLL

end

    
    % internal LLx (cross-validated likelihood, although was used for meta-params)
[LLs,~,~,fp] = fit0.eval_model(Robs, Xstim, XVi );
LLs(2) = fitS.eval_model(Robs, Xstim, XVi );
%LLs(3) = fit1.eval_model(Robs, Xstim, XVi );
LLs(3) = fit1S.eval_model(Robs, Xstim, XVi );
LLs(4) = fit2S.eval_model(Robs, Xstim, XVi );
LLs-fp.nullLL
    
    

