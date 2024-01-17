%% load a data set
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');

datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/qt0030/data000/data000');

datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

movie_path = '/Users/gdfield/Development/movie-xml2/BW-8-2-0.48-11111.xml';
datarun = load_java_movie(datarun, movie_path, datarun.triggers);
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params, 'keep_rf_coms', false);

%% set some params
corr_time = 0.01; % s


% hOS and OFF type 2
cell1 = 661;
cell2 = 1169;


% ON type 1 and type 2 Example 1
cell1 = 4;
cell2 = 6903;

% ON type 1 and type 2 Example 2
cell1 = 2461;
cell2 = 2821;

% Two ON type 1 cells: Example 1
cell1 = 4;
cell2 = 5567;

% Two ON type 2 cells: Example 1
cell1 = 6512;
cell2 = 6903;

% OFF type 1 and 2 cells
cell1 = 6391
cell2 = 6721;

% Two hOS cells: 10-31 datasets
cell1 = 6494;
cell2 = 6545;

cell1 = 6542;
cell2 = 7261;

cell_indices= get_cell_indices(datarun, [cell1 cell2]);


spikes_one = datarun.spikes{cell_indices(1)};
spikes_two = datarun.spikes{cell_indices(2)};


% compute ccf and plot
temp_options.dt = 0.005;
temp_options.offset = 0.2;
[ccf, ccf_time] = compute_ccf_fix(spikes_one, spikes_two, temp_options);
figure(1)
plot(ccf_time, ccf, 'k')
axis square



corr_spike_times = [];
for cnt = 1:length(spikes_one) 
    temp_ind = find(abs(spikes_two - spikes_one(cnt)) < corr_time);
    if ~isempty(temp_ind)
        %corr_spike_times = cat(1, corr_spike_times, spikes_two(temp_ind));
        corr_spike_times = cat(1, corr_spike_times, spikes_one(cnt));
    end
end

% bin up spikes for entire duration
spike_rate_ = histc(corr_spike_times, datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

% translate to spike times (with duplicates for multiple spikes per time bin)
spike_times_ = [];
for nn = 1:max(spike_rate_)
    spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
end




%%
frame_duration = datarun.stimulus.refresh_period * 0.001;

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
kept_frames = find(spike_times_ > 32);
spike_times_ = spike_times_(kept_frames);


T=text_waitbar;

stixel_num = height * width;
sta_depth = 32; % frames
clip_sum = zeros(stixel_num,sta_depth);
filt_frames = zeros(stixel_num,length(corr_spike_times));

% put spike times into a vector for this RGC instead of a cell array
temp_spike_times = spike_times_;
for fm = 1:length(temp_spike_times)
    
    T = text_waitbar(T,fm/length(temp_spike_times) - 0.01);   
    
    temp_mov = zeros(height,width,3,sta_depth);
    cntr = 0;
    for i= (temp_spike_times(fm) - sta_depth+1):temp_spike_times(fm)
        
       
        F = mvi.getFrame(i-1).getBuffer - 0.5;
        
        F = reshape(F,3,width,height);
        cntr = cntr + 1;
        temp_mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end

    % squeeze and reshape to sum the movie
    temp_clip = squeeze(temp_mov(:,:,1,:));
    clip_sum = clip_sum + reshape(temp_clip, [stixel_num, sta_depth]);

end

STA = clip_sum ./ length(temp_spike_times);
STA = reshape(STA, [40,80,sta_depth]);
STA = repmat(STA, [1,1,1,3]);
STA = permute(STA, [1,2,4,3]);

%% standard deconstruction of STA based on significant stixels
[sig_stixels, params, rf_strength_out] = significant_stixels(STA, 'thresh', 4);
temp_tc = time_course_from_sta(STA, sig_stixels);

temp_rf = reshape(squeeze(STA(:,:,1,:)), [stixel_num, sta_depth]) * temp_tc;
temp_rf = reshape(temp_rf, [40, 80, 3]);
norm_rf = norm_image(temp_rf);
imagesc(squeeze(norm_rf(:,:,1)))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
axis equal
axis tight
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

% filter the CORR RF
temp_params.radius = 1.0;
[filt_CORR_RF, params] = rf_filtered(norm_rf, temp_params);
imagesc(squeeze(filt_CORR_RF(:,:,1)))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
axis equal
axis tight
set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])


%% SVD based deconstruction

reshaped_STA = squeeze(STA(:,:,1,:));
reshaped_STA = reshape(reshaped_STA, [stixel_num, sta_depth]);

[U, S, V] = svd(reshaped_STA);

plot(diag(S), 'ko')
imagesc(reshape(U(:,1), [40,80]))

%% play STA

for fm = 1:size(STA, 4)
    imagesc(norm_image(squeeze(STA(:,:,:,fm))))
    axis equal
    pause
end

%%

% get sta summaries and filter
filt_params.radius = 1;
datarun = get_rfs_filtered(datarun, [cell1, cell2], 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');



%%
rf_one = datarun.stas.filt_rfs{cell_indices(1)};
rf_two = datarun.stas.filt_rfs{cell_indices(2)};
norm_rf_one = norm_image(rf_one);
norm_rf_two = norm_image(rf_two);

subplot(2,2,1)
imagesc(squeeze(norm_rf_one(:,:,1)));
colormap(brewermap([],'RdBu'))
caxis([0,1])
xlabel('cell A')
axis equal
axis tight

subplot(2,2,2)
imagesc(squeeze(norm_rf_two(:,:,1)));
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
xlabel('cell B')
axis equal
axis tight

subplot(2,2,3)
imagesc(squeeze(filt_CORR_RF(:,:,1)))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
xlabel('A x B')
axis equal
axis tight

subplot(2,2,4)
union_rfs = (norm_rf_one + norm_rf_two) ./2;
imagesc(squeeze(union_rfs(:,:,1)))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
xlabel('A + B')
axis equal
axis tight

print(1, '~/Desktop/corr_sta_comp.pdf', '-dpdf')


figure(2)


diff_rf = norm_image(squeeze(norm_rf(:,:,1)) - squeeze(union_rfs(:,:,1)));
imagesc(squeeze(diff_rf(:,:,1)))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
axis equal
axis tight

%% Compare time courses
cell1_tc = datarun.stas.time_courses{cell_indices(1)};
cell2_tc = datarun.stas.time_courses{cell_indices(2)};

figure(1); clf;
plot(temp_tc(3:32,1) ./ norm(temp_tc(3:32,1)), 'k')
hold on
plot(cell1_tc ./ norm(cell1_tc), 'r')
plot(cell2_tc ./ norm(cell2_tc), 'b')
hold off










