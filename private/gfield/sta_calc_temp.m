datapath = '/Users/gfield/Analysis/rat/2012-05-02-0/data001/data001';

datarun = load_data(datapath);
datarun = load_neurons(datarun);

movie_path = '/Users/gfield/Development/movie-xml2/BW-8-2-0.48-11111.xml';

% you need to define how many frames of the move to get
refresh_rate = 60; % this is the refresh rate of the video display in Hz.
num_frames = floor(datarun.duration * refresh_rate);
frame_bins = datarun.triggers(2:end):(1/refresh_rate):datarun.duration;

temp_index = get_cell_indices(datarun, 18);
spike_times = datarun.spikes{temp_index};

num_frames = 50000;

[mov,height, width, duration, refresh] = get_movie(movie_path, datarun.triggers, num_frames);

mov = squeeze(mov(:,:,1,:));
mov = permute(mov, [3, 1, 2]);
mov = reshape(mov, [num_frames, 40*80]);





rgc_STA = zeros(height, width, STA_length);
for Vtrig = 3:length(datarun.triggers)-1 
    [mov,height,width,duration,refresh] = get_movie(movie_path, datarun.triggers(Vtrig-1:Vtrig), 200);
    mov = squeeze(mov(:,:,1,:));
    mov = 2*mov - 1;
    
    temp_times = spike_times(spike_times >= datarun.triggers(Vtrig) & spike_times < datarun.triggers(Vtrig +1));
%    temp_bins = datarun.triggers(Vtrig):(1/refresh_rate):datarun.triggers(Vtrig +1);
    
    binned_spike_times = histcounts(temp_times,100);
    spike_triggers = find(binned_spike_times > 0);
    
    for fm = 1:length(spike_triggers)
        temp_mov = mov(:,:,100+spike_triggers(fm)-STA_length+1:100+spike_triggers(fm));
        temp_mov = temp_mov * binned_spike_times(spike_triggers(fm));
        rgc_STA = rgc_STA + temp_mov;
    end
end

%%
[mvi] = load_movie(movie_path, datarun.triggers);

F = mvi.getFrame(0).getBuffer;

F = reshape(F,3,width,height);
mov(:,:,:,i) = permute(F,[3 2 1]);



%%



% size the movie to check its size.
spike_times = datarun.spikes{1};
spike_times = spike_times(spike_times > datarun.triggers(2));


% bin spike times
bin_spike_times = histc(spike_times, frame_bins);

spike_triggers = find(bin_spike_times > 0);

STA_length = 24;


rgc_STA = zeros(height, width, STA_length);
for fm = 1:length(spike_triggers)
    if spike_triggers(fm) < STA_length
        continue
    else
        temp_mov = mov(:,:,spike_triggers(fm)-STA_length+1:spike_triggers(fm));
        temp_mov = temp_mov * bin_spike_times(spike_triggers(fm));
    end
    rgc_STA = rgc_STA + temp_mov;
end

rgc_STA = rgc_STA ./ length(spike_times);

for fm = 1:STA_length
    image(norm_image(squeeze(rgc_STA(:,:,fm))))
    pause(0.2)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
movie_path = '/Users/gfield/Development/movie-xml2/BW-15-1-0.48-11111-52x40-60.35.xml';
datarun = load_java_movie(datarun, movie_path, datarun.triggers);

frame_bins = datarun.triggers(1:end):(1/refresh_rate):datarun.duration;

refresh_rate = 60.35; % this is the refresh rate of the video display in Hz.
num_frames = floor(datarun.duration * refresh_rate);

[gmov,height,width,duration,refresh] = get_movie_from_datarun(datarun, num_frames);
%gmov = squeeze(gmov(:,:,1,:));
gmov = 2*gmov - 1;

% size the movie to check its size.
spike_times = datarun.spikes{1};


% bin spike times
bin_spike_times = histcounts(spike_times, frame_bins);

spike_triggers = find(bin_spike_times > 0);

STA_length = 24;
rgc_STA = zeros(height, width, STA_length);
for fm = 1:length(spike_triggers)
    if spike_triggers(fm) < STA_length
        continue
    else
        temp_mov = gmov(:,:,spike_triggers(fm)-STA_length+1:spike_triggers(fm));
        temp_mov = temp_mov * bin_spike_times(spike_triggers(fm));
    end
    rgc_STA = rgc_STA + temp_mov;
end

rgc_STA = rgc_STA ./ length(spike_times);

for fm = 1:STA_length
    image(norm_image(squeeze(rgc_STA(:,:,fm))))
    pause(0.2)
end

