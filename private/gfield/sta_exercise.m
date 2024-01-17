% This data should be on your local drive, so you can change the path on
% the line below to point it your local version: this will save download
% time.
data_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-11-27-0/data001/data001';

% load all the relevant stuff
datarun = load_data(data_path);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);

% find the cell ID for the first ON type1 cell.
info(datarun)

% Plot the receptive field (RF) of this cell.
temp_cell_ID = 271;
plot_rf(datarun, temp_cell_ID)

% Now you are going to calculate the spike-triggered average (STA) of 
% this cell

% first get the spike times.
temp_index = get_cell_indices(datarun, temp_cell_ID); % gets the index for cell 271
temp_spike_times = datarun.spikes{temp_index}; % gets the spike times for this cell

% To compute the STA, you will also need the stimulus that was delivered.
% Here is how to get it.
% First define the path to the xml file that sets has white noise movie parameters.
% Note: This xml sets things like the size of the white noise checkers and how
% often they refresh
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-10-1-0.48-11111-60x60-60.35.xml';

% you need to define how many frames of the move to get
refresh_rate = 60.35; % this is the refresh rate of the video display in Hz.
num_frames = floor(datarun.duration * refresh_rate);

% Note the number of frames is 337,175 and each frame of the movie is
% described by a 60x60 grid of pixels with 3 values at each pixel to define
% the color (red, green, and blue intensities for the display primaries).
% Thus, the movie is an array that is 337,175 x 60 x 60 x 3. This is too
% large to load in the Matlab's memory all at once. So we are just going to
%  load in the first 100,000 frames (note, you might have to change the max array
%  size in matlab to handle even this).

% You'll need to connect to the server with VPN (dusom_fieldlab) to get to this.
num_load_frames = 100000;
[mov,height,width,duration,refresh] = get_movie(movie_path, datarun.triggers, num_load_frames);

% size the movie to check its size.
size(mov)

% The movie is black and white (BW), meaning that each pixel can be only
% black or white, so we can make this matrix a little smaller, by reducing
% the color dimension
mov = squeeze(mov(:,:,1,:));

% now size the movie:
size(mov)

% OK, so now you have the visual stimulus (contained in 'mov') and the spike
% times for one of the neurons (contained in 'spike_times'). Your exercise is
% to write code to compute the spike triggered average and the spike
% triggered covariance for this neuron. See the Schwartz2006 paper for how to do each of
% these.  Let me know if you get stuck or questions come up along the way.

%% Bin the spikes to correspond with the movie

% define a vector of frame times
frame_times = 0:((1/refresh_rate) * 2):((1/refresh_rate) * 2)*num_load_frames;

% bin the spikes according to the frame times (this is the number of spikes
% happening on each frame)
[spike_counts, ~] = histcounts(temp_spike_times, frame_times);

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









