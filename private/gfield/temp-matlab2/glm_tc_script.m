%% Load data

datapath = '/Users/gfield/Analysis/2017-01-16-0/data006_KR/data006_KR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';
%moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35.xml';

%datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-06-24-0/data005_MR/data005_MR';
%moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
datapath = '/Volumes/Disk2/Data/2019-06-24-0/data005_MR/data005_MR';
moviepath = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-53x40-60.35.xml';


datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

cell_types = {'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
cell_indices = get_cell_indices(datarun, cell_types);
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, cell_types, 'marks_params', marks_params);

% get movie
frame_num = 1600;
[mov,height,width,duration,refresh] = get_movie(moviepath, datarun.triggers, frame_num);


%% Choose neuron, and plot raster


for rgc = 1:length(cell_indices)

    temp_spikes = datarun.spikes{cell_indices(rgc)};

    % get center of mass and of STA and cut out an ROI to reduce dim of stim
    temp_com = abs([0 40] - round(datarun.vision.sta_fits{cell_indices(rgc)}.mean)) +[1 1];
    temp_size = 5*sqrt(prod(datarun.vision.sta_fits{cell_indices(rgc)}.sd));
    window_size = round(temp_size);
    bg_height = floor(temp_com(2) - window_size);
    end_height = floor(temp_com(2) + window_size);
    bg_width = floor(temp_com(1) - window_size);
    end_width = floor(temp_com(1) + window_size);
    wind_dim = window_size*2 + 1;

    new_mov = squeeze(mov(bg_height:end_height, bg_width:end_width,1,:));
    new_mov = reshape(new_mov, [(wind_dim * wind_dim), frame_num]);
    new_mov = permute(new_mov, [2 1]) -0.5;
    dt = refresh./ 1000;

    frame_time = refresh * frame_num /1000;
    %cut_spikes = temp_spikes(temp_spikes > datarun.triggers(1) &  temp_spikes < datarun.triggers(frame_num/100));

    %% Fit GLM model
    NIM_params.nLags = 30; %number of time lags for estimating stimulus filters
    NIM_params.dt = refresh./ 1000;
    NIM_params.up_samp_fac = 8; % temporal up-sampling factor applied to stimulus 
    NIM_params.tent_basis_spacing = 8; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
    NIM_params.silent = 0;
    NIM_params.binsize = NIM_params.dt/NIM_params.up_samp_fac;
    NIM_params.model_to_fit = [1];
    NIM_params.wind_dim = wind_dim;

    NIM_params.NX = wind_dim;
    NIM_params.NY = wind_dim;

    temp_triggers = datarun.triggers(1:(frame_num/100)+1);
    cut_spikes = temp_spikes(temp_spikes > datarun.triggers(1) &  temp_spikes < datarun.triggers(frame_num/100));

    fit_struct = NIM_wrapper_glm_only(new_mov, cut_spikes, temp_triggers, NIM_params);
    ndf0_fits{rgc} = fit_struct;
    
    
    % extract TC and RF from GLM fit w/ spike feedback
    glm_filt = fit_struct.fitS.subunits.filtK;
    glm_sta =  reshape(glm_filt, [NIM_params.nLags, NIM_params.NX*NIM_params.NY]);
    glm_sta = permute(glm_sta, [2, 3, 1]);
    glm_sta = repmat(glm_sta, [1 1 1 3]);
    glm_sta = permute(glm_sta, [1 2 4 3]);
    glm_filt = reshape(glm_filt, [NIM_params.nLags, NIM_params.NX*NIM_params.NY]);
    glm_std = std(glm_filt(:));
    [glm_i, glm_j] = find(abs(glm_filt) > 2.5*glm_std);
    marks_matrix = zeros(size(glm_filt));
    marks_matrix(glm_i, glm_j) = 1;
    marks_matrix = max(marks_matrix);
    marks_matrix = reshape(marks_matrix, [NIM_params.NX, NIM_params.NY]);
    sig_stix_glm = logical(marks_matrix);
    temp_tc_glm = time_course_from_sta(glm_sta, sig_stix_glm);
    temp_rf_glm = rf_from_sta(glm_sta, 'sig_stixels', sig_stix_glm);
    

    % extract TC and RF from ~LNP (GLM fit w/o spike feedback)
    lnp_filt = fit_struct.fit0.subunits.filtK;
    lnp_sta = reshape(lnp_filt, [NIM_params.nLags, NIM_params.NX, NIM_params.NY]);
    lnp_sta = permute(lnp_sta, [2, 3, 1]);
    lnp_sta = repmat(lnp_sta, [1 1 1 3]);
    lnp_sta = permute(lnp_sta, [1 2 4 3]);
    lnp_filt = reshape(lnp_filt, [NIM_params.nLags, NIM_params.NX*NIM_params.NY]);
    lnp_std = std(lnp_filt(:));
    [lnp_i, lnp_j] = find(abs(lnp_filt) > 2.5*lnp_std);
    marks_matrix = zeros(size(lnp_filt));
    marks_matrix(lnp_i, lnp_j) = 1;
    marks_matrix = max(marks_matrix);
    marks_matrix = reshape(marks_matrix, [NIM_params.NX, NIM_params.NY]);
    sig_stix_lnp = logical(marks_matrix);     
    temp_tc_lnp = time_course_from_sta(lnp_sta, sig_stix_lnp);
    temp_rf_lnp = rf_from_sta(lnp_sta, 'sig_stixels', sig_stix_lnp);

    
    figure(1); clf;  
    rev_indices = NIM_params.nLags:-1:1;
    plot(temp_tc_glm(rev_indices,1) ./ ext(temp_tc_glm(rev_indices,1)), 'k')
    hold on
    plot(temp_tc_lnp(rev_indices,1) ./ ext(temp_tc_glm(rev_indices,1)), 'r')
    plot(datarun.stas.time_courses{cell_indices(rgc)} ./ ext(datarun.stas.time_courses{cell_indices(rgc)}), 'b');
    hold off
    title(num2str(datarun.cell_ids(cell_indices(rgc))));
    drawnow
    
    %fit_struct.fit0.display_model
    %fit_struct.fit2S.display_model

end





figure(2); clf;
semilogy(diag(lnp_S).^2 ./ sum(diag(lnp_S).^2), 'ro')
hold on
semilogy(diag(glm_S).^2 ./ sum(diag(glm_S).^2), 'ko')
hold off




