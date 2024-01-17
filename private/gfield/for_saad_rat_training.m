% 2012-05-02-0 data001
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-05-02-0/data001/data001';
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-8-2-0.48-11111.xml';

% 2012-10-15-0 data000
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/data000-map/data000-map';
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-8-2-0.48-11111.xml';

% 2012-10-31-0 data000
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-31-0/data000-1800-7200/data000-map/data000-map';
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-8-2-0.48-11111.xml';



datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

%num_frames = 200000;
num_frames = datarun.duration * 60;
[temp_mov,~,~,dur,frame_refresh] = get_movie_binary(movie_path, datarun.triggers, num_frames); %num of unique frames
temp_mov = squeeze(temp_mov(:,:,1,:)) * 2 - 1.0;

temp_bins = [datarun.triggers(1):frame_refresh:(datarun.triggers(1) + (num_frames * frame_refresh))] ./1000;

temp_spikes = datarun.spikes{2};
binned_spikes = histcounts(temp_spikes, temp_bins);

nLags = 24;
temp_sta = zeros(40,80,nLags);

for fm = nLags:num_frames    
    temp_sta = temp_sta + double((binned_spikes(fm) * temp_mov(:,:, (fm - nLags+1):fm)));
end

for fm = 1:nLags
    imagesc(squeeze(temp_sta(:,:,fm)))
    axis square
    pause(0.25)
end

cd ~/Desktop
% save elements for Saad
model_structure.spike = datarun.spikes;
model_structure.triggers = datarun.triggers;
model_structure.cell_ids = datarun.cell_ids;
model_structure.cell_types = datarun.cell_types;
model_structure.movie = temp_mov;

save 2012-10-31-0-data000 model_structure


%%

load 2012-10-31-0-data000

num_frames = size(model_structure.movie, 3);
frame_refresh = (model_structure.triggers(2) - model_structure.triggers(1)) ./ 50;
temp_bins = [model_structure.triggers(1):frame_refresh:(model_structure.triggers(1) + (num_frames * frame_refresh))];

% look at the RF for the a particular cell
cell_index = 2;
temp_spikes = model_structure.spike{cell_index};
binned_spikes = histcounts(temp_spikes, temp_bins);

% number of frames back to compute the STA
nLags = 24;
% prealocate the STA
temp_sta = zeros(40,80,nLags);

for fm = nLags:num_frames    
    temp_sta = temp_sta + (double((binned_spikes(fm) * model_structure.movie(:,:, (fm - nLags+1):fm))));
end

% look at the 24 frames of the STA
for fm = 1:nLags
    imagesc(squeeze(temp_sta(:,:,fm)))
    axis square
    pause(0.25)
end


