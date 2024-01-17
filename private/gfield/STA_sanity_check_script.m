cd ~/Dropbox/shares/NIM_work/NIM-tests/cell5327/
load spikes
load mov


%% VECTORIZED STA calculation
nLags = 24;  
dt = 16.5975 ./ 1000;
binsize = dt;

% Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
[NFRAMES, nPix] = size(mov);
NX = 11;  % can handle 2-d stim, just need to specify in stim params:

Xstim = zeros(size(mov,1), size(mov,2) * nLags);
for fm = 1:nLags
    if fm == 1
        Xstim(:,fm-1+(1:nLags:(nPix*nLags))) = mov;
    else
        Xstim(:,fm-1+(1:nLags:(nPix*nLags))) = [zeros(fm, nPix); mov(1:end-fm,:)];           
    end
end
         
bin_edges = 0:dt:NFRAMES*dt;
binned_spikes = histcounts(spikes,bin_edges);
spike_triggers = find(binned_spikes > 0);
STE = Xstim(spike_triggers, :);
STA = mean(STE);
STA = reshape(STA, [nLags, NX, NX]);

for fm = 1:nLags
    imagesc(squeeze(STA(fm,:,:)))
    pause(0.2)
end


%% LOOP STA calculation
nLags = 24;  
dt = 16.5975 ./ 1000;
binsize = dt;

bin_edges = 0:dt:NFRAMES*dt;
binned_spikes = histcounts(spikes,bin_edges);

spike_triggers = find(binned_spikes > 0);
spike_triggers = spike_triggers(spike_triggers > nLags);


STA = zeros(nLags,size(mov, 2)); 
for fm = 1:length(spike_triggers)
    tmp_movie_clip = mov(spike_triggers(fm)-nLags+1:spike_triggers(fm),:);
    STA = STA + tmp_movie_clip * binned_spikes(spike_triggers(fm));
end

STA = STA ./ length(spikes);

STA = reshape(STA, [nLags, NX, NX]);

for fm = 1:nLags
    imagesc(squeeze(STA(fm,:,:)))
    pause(0.2)
end

%%
