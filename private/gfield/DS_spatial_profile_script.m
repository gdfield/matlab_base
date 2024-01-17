% load data, compute RFs and TCs
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

filt_params.radius = 1;

datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

ds_right = [257 1097 1683 1895 2898 3512 3736 4069 4157 4846 5632 6321 6751 6797];
ds_left = [301 2449 3934 4173 4789 5719 6332 7308];
ds_up = [708 3815 6229];
ds_down = [333 438 467 1595 1685 2042 3636 3842 4353 4731 4985 5150 5702];

offt3 = [62,991,1156,4234,4278,4487,5733,6286,6931];
offt5 = [1246,2253,3695,5116,6260];

%% load 2012-10-31-0

datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

filt_params.radius = 1;

datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

%% load 2012-10-10-1

datarun = load_data('/Users/gdfield/Analysis/2012-10-10-1/data001-3600-7200s/data001-map/data001-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

filt_params.radius = 1;

datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

%%
% book keeping
cell_type_num = get_cell_type_nums(datarun, 'off type4');
cell_indices = get_cell_indices(datarun, {cell_type_num});

num_rgcs = length(cell_indices);
singval_ratios = zeros(1,num_rgcs);
for rgc = 1:num_rgcs
    
    cell_id = datarun.cell_ids(cell_indices(rgc))
    
    % get the space-time STA
    temp_sta = datarun.stas.stas{cell_indices(rgc)};
    temp_sta = squeeze(temp_sta);
    
    % apply a spatial filter to STA
    if 1
        filt_sta = squeeze(zeros(size(datarun.stas.stas{1})));
        temp_params.radius = 1;
        for fm = 1:datarun.stas.depth
            temp_fm = rf_filtered(temp_sta(:,:,fm), temp_params);
            filt_sta(:,:, fm) = temp_fm;
        end
    end
    filt_sta = reshape(filt_sta, [3200, 30]);
    
    % squeeze and reshape the sta for SVD
    temp_sta = reshape(filt_sta, [3200,30]);

    % perform SVD
    [U, S, V] = svd(filt_sta(:,11:30));

    % calculate and store singular value ratios 
    singvals = diag(S);
    si = singvals.^2 ./ sum(singvals.^2);
    
    
    model_sta = U(:,1)*S(1)*V(:,1)';
    y = reshape(filt_sta(:,11:30), 1,[]);
    x = reshape(model_sta, 1, []);
    p = polyfit(x, y, 1);
    yfit = polyval(p,x);
    yfit =  p(1) * x + p(2);   
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal
    
    
%     % compute a "noise" sta and do SVD on that.
%     noise_mean = median(temp_sta(:));
%     noise_std = std(temp_sta(:));
%     noise_sta = normrnd(noise_mean, noise_std, size(temp_sta));
%     [U_noise, S_noise, V_noise] = svd(noise_sta);

    % plot stuff
    num_dims = 4;
    spaces = zeros(40, 80, num_dims);
    temporals = zeros(30, num_dims);
    score_thresh = 3.5;
    for dims = 1:num_dims
        spaces(:, :, dims) = reshape(U(:,dims), [40, 80]);
        figure(1)
        subplot(2, 5, dims+1)
        temp_im = (norm_image(spaces(:,:, dims)));
        imagesc(matrix_scaled_up(squeeze(temp_im(:,:,1)),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        axis off
        axis square

        %temporals(:,dims) = V(:,dims);
        %figure(2)
        subplot(2,5,dims+6)
        plot(V(:,dims), 'k')
        axis('tight')
        axis square
                
        subplot(2,5,6)
        plot(si, 'k.', 'MarkerSize', 8)
        axis([0 16 0 1])
        axis square
    
    end
        
    %figure(4)
    subplot(2,5,1)
    temp_rf = datarun.stas.filt_rfs{cell_indices(rgc)};
    norm_rf = norm_image(temp_rf);
    new_rf = squeeze(norm_rf(:,:,1));
    imagesc(matrix_scaled_up(new_rf,8))
    axis off
    axis square
    
    colormap(brewermap([],'RdBu'))
    caxis([0,1])
    
    si_1 = si(1)

    
    drawnow
    pause
end

%%
print(1, '/Users/gdfield/Desktop/test-fig.pdf', '-dpdf')

%%
temp_bins = 0:0.05:1;
on1_hist = hist(on1_singval_ratios, temp_bins);
on2_hist = hist(on2_singval_ratios, temp_bins);

plot(temp_bins, on1_hist, 'r', temp_bins, on2_hist, 'k')

figure(2)
bar(temp_bins, on1_hist, 'r')
hold on
bar(temp_bins, on2_hist, 'k')
hold off
print(2, '/Users/gdfield/Desktop/hist-fig.pdf', '-dpdf')


%% compute a spatial STC after filtering the stimulus in time
% load java object of movie
datarun = load_java_movie(datarun); 
    
% refresh times
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;


temp_cell_id = 2461;
temp_cell_index = get_cell_indices(datarun, temp_cell_id);
temp_sta = datarun.stas.stas{temp_cell_index};

 % squeeze and reshape the sta for SVD
temp_sta = squeeze(temp_sta);
temp_sta = reshape(temp_sta, [3200,30]);

% perform SVD
[U, S, V] = svd(temp_sta);

% get rank-1 time estimate
tc_estimate = V(:,2);
%tc_estimate = tc_estimate - mean(tc_estimate);
plot(tc_estimate)
sta_tc_length = length(tc_estimate);


%% filter the movie by the TC

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
mvi = load_movie('/Users/gdfield/Development/movie-xml2/BW-8-2-0.48-11111.xml', datarun.triggers);
% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

total_frame_number = datarun.duration ./ (datarun.stimulus.refresh_period * 0.001);

% only keep spikes that are more than 30 frames from the beginning of the
% stimulus
first_kept_frame = find(spike_times{1} > 30, 1);
spike_times{1} = spike_times{1}(first_kept_frame:end);


T=text_waitbar;

clip_sum = zeros(3200,30);
filt_frames = zeros(3200,length(spike_times));

% put spike times into a vector for this RGC instead of a cell array
temp_spike_times = spike_times{1};
for fm = 1:length(temp_spike_times)
    
    T = text_waitbar(T,fm/length(temp_spike_times) - 0.01);   
    
    temp_mov = zeros(height,width,3,sta_tc_length);
    cntr = 0;
    for i= (temp_spike_times(fm) - sta_tc_length+1):temp_spike_times(fm)
        
        F = mvi.getFrame(i-1).getBuffer;
        
        F = reshape(F,3,width,height);
        cntr = cntr + 1;
        temp_mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end

    % squeeze and reshape to sum the movie
    temp_clip = squeeze(temp_mov(:,:,1,:));
    clip_sum = clip_sum + reshape(temp_clip, [3200, sta_tc_length]);
    
    
%     % filter the frame
     temp_filt_frames = reshape(temp_mov(:,:,1,:), [3200,sta_tc_length]) * tc_estimate;
 
     filt_frames(:,fm) = temp_filt_frames;

end

STA = reshape(clip_sum, [40,80,30]);

for fm = 1:length(STA)
    imagesc(STA(:,:,fm))
    pause
end

% look at the average of frame_keeper
st_mean = mean(frame_keeper,2);
st_mean = reshape(st_mean, [40, 80]);
imagesc(st_mean); colormap gray

% look at the covariance of frame_keeper
stv = var(filt_frames, [],2);

stv = reshape(stv, [40,80]);
imagesc(stv); colormap jet















