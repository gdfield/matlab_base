%% compute spike-trigger covariance of OS cells

datapath = '/Volumes/gdf/rat-data/2012-10-15-0/YASS/data000/data000';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'OFF vOS');

marks_params.thresh = 3.5;
datarun = get_sta_summaries(datarun, 'OFF vOS', 'marks_params', marks_params);
vOS_indices = get_cell_indices(datarun, 'OFF vOS');

%%

movie_path = '/Users/gfield/Library/CloudStorage/Dropbox/Work/Documents/Development/movie-xml2/BW-8-2-0.48-11111.xml';
refresh_rate = 60; % this is the refresh rate of the video display in Hz.
num_frames = floor(datarun.duration * refresh_rate);
frame_bins = datarun.triggers(1):(1/refresh_rate):datarun.duration;

spike_times = datarun.spikes{vOS_indices(1)};

num_frames = 50000;

[mov,height, width, duration, refresh] = get_movie(movie_path, datarun.triggers, num_frames);

mov = squeeze(mov(:,:,1,:));
mov = permute(mov, [3, 1, 2]);
mov = reshape(mov, [num_frames, 40*80]);


%%

movie_path = '/Users/gfield/Library/CloudStorage/Dropbox/Work/Documents/Development/movie-xml2/BW-8-2-0.48-11111.xml';
spike_times = datarun.spikes{vOS_indices(1)};

datarun = load_java_movie(datarun, movie_path, datarun.triggers);

refresh_rate = 120; % this is the refresh rate of the video display in Hz.
num_frames = floor(datarun.duration * refresh_rate);

num_load_frames = 100000;
[mov,height,width,duration,refresh] = get_movie(movie_path, datarun.triggers, num_load_frames);

mov = squeeze(mov(:,:,1,:));
size(mov)

mov_interval = datarun.stimulus.interval;
frame_times = datarun.triggers(1):((1/refresh_rate) * mov_interval):((1/refresh_rate) * mov_interval)*num_load_frames;

[spike_counts, ~] = histcounts(spike_times, frame_times);

% define the length of the sta
sta_length = 24;

% identify frames that have alteast 1 spike;
triggered_frames = find(spike_counts > 0);

% only keep the frames that occur after sta_length
triggered_frames = triggered_frames(triggered_frames >= sta_length);

temp_STA = zeros(height, width, sta_length);
for fm = 1:length(triggered_frames)
    begin_frame = triggered_frames(fm) - sta_length + 1;
    end_frame = triggered_frames(fm);
    temp_STA = temp_STA + (mov(:,:,begin_frame:end_frame) * spike_counts(triggered_frames(fm))); % need to multiply by the number of spikes for the STA;
end
temp_STA = temp_STA ./ sum(spike_counts(triggered_frames));


% view the STA
for fm = 1:sta_length
    imagesc(temp_STA(:,:,fm))
    pause(0.25)
    colormap gray
end


%%

datarun = get_sta_summaries(datarun, 'OFF type1', 'marks_params', marks_params);
temp_indices = get_cell_indices(datarun, 'OFF type1');

temp_spike_times = datarun.spikes{temp_indices(1)};

mov_interval = datarun.stimulus.interval;
refresh_rate = 120 / mov_interval;
%[gmov,height,width,duration,refresh] = get_movie_from_datarun(datarun, num_frames);
[gmov,height, width, duration, refresh] = get_movie(movie_path, datarun.triggers, 432300-1);
%gmov = squeeze(gmov(:,:,:,:));
gmov = 2*gmov - 1;

samples_per_interval = 50;
uniform_times = uniform_subsample_triggers(datarun.triggers, samples_per_interval);

hist_spikes = histcounts(temp_spike_times, uniform_times);

frames_back = 24;
[STA, STV] = compute_sta_stv_binned(gmov, hist_spikes, frames_back);

min_sta = min(STA(:));
max_sta = max(STA(:));
% check if abs(min) is larger than max
if abs(min_sta) < max_sta
    ext_sta = abs(min_sta);
else
    ext_sta = max_sta;
end
norm_STA = STA ./ ext_sta;
%make the STA consistent with image format
norm_STA = (norm_sta +1 ) ./2;


for fm = 1:frames_back
    image(norm_STA(:,:,:,fm))
    colormap gray
    pause(0.2)
end

reshaped_STA = reshape(squeeze(norm_STA(:,:,1,:)), [3200, frames_back]);
[spatial_kerns, eigen_vals, temporal_kerns] = svd(reshaped_STA);







% create a 

num_triggers = length(datarun.triggers);
frame_bins = zeros(duration+1, 1);
frame_bins(1)= datarun.triggers(1);
for tg = 1:num_triggers-1
    temp_bins = datarun.triggers(tg):refresh/1000:datarun.triggers(tg+1);
    size(temp_bins)
    frame_bg = 1+((tg-1) * 50);
    frame_end = tg * 50;
    frame_bins(frame_bg+1:frame_end+1) = temp_bins(2:51);
end




