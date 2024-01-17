% define data path and movie
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-07-18-0/data009-control/data009-control';
moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-52x40-60.35_xoffset.xml';


%% load data

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');


type_num = get_cell_type_nums(datarun, 'OFF brisk transient');
cell_ids = datarun.cell_types{type_num}.cell_ids;
cell_indices = get_cell_indices(datarun, cell_ids);

%% this where you would define how much data to use in terms of the number of triggers in the movie.
%num_trigs_for_calc = 700;
num_trigs_for_calc = length(datarun.triggers)-1;

% bin the time according to the trigger times.
bin_edges = datarun.triggers(1);
for trg = 1:num_trigs_for_calc
    diff_trig_time = datarun.triggers(trg+1) - datarun.triggers(trg);
    est_frame_length = diff_trig_time./100;
    tmp_times = datarun.triggers(trg)+est_frame_length : est_frame_length : datarun.triggers(trg+1);    
    bin_edges = cat(2, bin_edges, tmp_times);
end

%% Choose neuron
cell_id = datarun.cell_ids(cell_indices(2)); % example cell
datarun = get_sta_summaries(datarun, cell_id);
cell_index = get_cell_indices(datarun, cell_id);
temp_spikes = datarun.spikes{cell_index};

bin_spikes = histcounts(temp_spikes, bin_edges);

% load in the movie from the move field
num_frames = num_trigs_for_calc * 100; %100 frames per trigger
[mov,height,width,duration,refresh] = get_movie(moviepath, datarun.triggers, num_frames);
mov = squeeze(mov(:,:,1,:));
mov = reshape(mov, [(height * width), num_frames]) - 0.5;


% calculate the STA in tmp_sta
num_lags = 30;
tmp_sta = zeros(height * width, num_lags);
shift_bin_spikes = bin_spikes;
last_bin = length(shift_bin_spikes);
for nlag = 1:num_lags
    singles = shift_bin_spikes == 1;
    doubles = shift_bin_spikes == 2;
    tripples = shift_bin_spikes == 3;
    quads = shift_bin_spikes == 4;
    quints = shift_bin_spikes == 5;

    temp_frame = sum(mov(:,singles'), 2) + 2*sum(mov(:,doubles), 2) + 3*sum(mov(:,tripples), 2) + 4*sum(mov(:,quads), 2) + 5*sum(mov(:,quints), 2);
    tmp_sta(:, num_lags - nlag + 1) = temp_frame;
    
    % shift binned spikes to earlier time frames to collect frames at
    % relevant lags
    shift_bin_spikes = circshift(shift_bin_spikes, -nlag);
    shift_bin_spikes(last_bin) = 0;
end

tmp_sta = tmp_sta ./ sum(bin_spikes);

[U, S, V] = svd(tmp_sta); % get spatial and temporal RFs using Singular Value Decomposition

svd_dim = 2; % you may have to make this 1 instead of 2 depending on the cell you choose

figure(1)
imagesc(reshape(U(:,svd_dim), height, width)); % RF computed from limited data

figure(2); clf;
plot(V(:,svd_dim) ./ ext(V(:,svd_dim))) % TC computed from limited data



