%%
datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/data004/data004');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

filt_params.radius = 1;

datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% DS cell
temp_cell_id = 361;

% high SNR comparison
%temp_cell_id = 1052;
temp_cell_index = get_cell_indices(datarun, temp_cell_id);



% load java object of movie
datarun = load_java_movie(datarun); 
    
% refresh times
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;


% bin up spikes for entire duration
spike_rate_ = histc(datarun.spikes{temp_cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

% translate to spike times (with duplicates for multiple spikes per time bin)
spike_times_ = [];
for nn = 1:max(spike_rate_)
    spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
end

% put into storage variables
spike_rate(1,:) = spike_rate_';
spike_times{1} = sort(spike_times_);


%%

frame_duration = datarun.stimulus.refresh_period * 0.001;

% initialize frame_keeper to zero
frame_keeper = zeros(3200, length(spike_times{1}));
cntr1 = 0;

% get movie object
mvi = load_movie('/Users/gdfield/Development/movie-xml2/BW-gaussian-8-4-0.16-11111.xml', datarun.triggers);

mvi = load_movie('/Users/gdfield/Development/movie-xml2/BW-8-4-0.48-11111.xml', datarun.triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

total_frame_number = datarun.duration ./ (datarun.stimulus.refresh_period * 0.001);

% only keep spikes that are more than 30 frames from the beginning of the
% stimulus
first_kept_frame = find(spike_times{1} > 16, 1);
spike_times{1} = spike_times{1}(first_kept_frame:end);


for fm = 1:(datarun.duration * datarun.stimulus.monitor_refresh ./ datarun.stimulus.interval);
    F = (mvi.getFrame(fm-1).getBuffer - 0.5) * 2;
    
     % running tally of the variance
    if fm == 1
        Mk = F;
        Sk = 0;
    else
        Mkm1 = Mk;
        Mk = Mk + ((F - Mk)/fm);
        Sk = Sk + ((F - Mkm1).*(F - Mk));
    end
end


stim_ensemble_var = Sk/(fm - 1); 

stim_ensemble_var = reshape(stim_ensemble_var,3,width,height);
stim_ensemble_var = permute(stim_ensemble_var,[3 2 1]);
stim_ensemble_var = stim_ensemble_var(:,:,1);

%%

T=text_waitbar;

clip_sum = zeros(3200,16);
filt_frames = zeros(3200,length(spike_times));

% put spike times into a vector for this RGC instead of a cell array
temp_spike_times = spike_times{1};
ste = zeros((80*40*16), length(temp_spike_times));
for fm = 1:length(temp_spike_times)
    
    T = text_waitbar(T,fm/length(temp_spike_times) - 0.01);   
    
    temp_mov = zeros(height,width,3,16);
    cntr = 0;
    for i= (temp_spike_times(fm) - 16+1):temp_spike_times(fm)
        
        F = (mvi.getFrame(i-1).getBuffer - 0.5) * 2;
        
        F = reshape(F,3,width,height);
        cntr = cntr + 1;
        temp_mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end

    % squeeze and reshape to sum the movie
    temp_clip = squeeze(temp_mov(:,:,1,:));
    temp_clip = reshape(temp_clip, [], 1);
    
    % running tally of the variance
    if fm == 1
        Mk = temp_clip;
        Sk = 0;
    else
        Mkm1 = Mk;
        Mk = Mk + ((temp_clip - Mk)/fm);
        Sk = Sk + ((temp_clip - Mkm1).*(temp_clip - Mk));
    end

    ste(:,fm) = temp_clip;
    
end

%%
% stv = var(ste, 0, 2);
% stv = stv - 0.0256;
% stv = reshape(stv, [40, 80, 16]);

%%
sta = mean(ste,2);
sta = reshape(sta, [40, 80, 16]);

%%
stv = Sk/(fm - 1); 
stv = reshape(stv, [40, 80, 16]);
stv = stv - repmat(stim_ensemble_var, [1,1,16]);



%% inspect STA
for fm = 1:16
    imagesc(sta(:,:,fm))
    %colormap(brewermap([],'RdBu'))       
    %caxis([0,1])        
    pause
end



%% inspect STV
for fm = 1:16
    imagesc(stv(:,:,fm))
    %colormap(brewermap([],'RdBu'))       
    %caxis([0,1])        
    pause
end

%%

STA = reshape(clip_sum, [40,80,16]);

% look at the average of frame_keeper
st_mean = mean(frame_keeper,2);
st_mean = reshape(st_mean, [40, 80]);
imagesc(st_mean); colormap gray

% look at the covariance of frame_keeper
stv = var(filt_frames, [],2);

stv = reshape(stv, [40,80]);
imagesc(stv); colormap jet




