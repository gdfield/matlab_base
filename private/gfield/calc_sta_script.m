%% STA calculation in matlab
data_path = '/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map';
movie_path = '/Users/gdfield/Development/movie-xml2/BW-8-2-0.48-11111.xml';

% specify cell
temp_cell_index = 2;

%specify the number of frames for the STA.
sta_frame_num = 15;

%Get datarun and java movie
datarun = load_data(data_path);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_java_movie(datarun); 
    
% refresh times and frame duration
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
frame_duration = datarun.stimulus.refresh_period * 0.001;

% get movie object
mvi = load_movie(movie_path, datarun.triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);
total_frame_number = datarun.duration ./ (datarun.stimulus.refresh_period * 0.001);

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

% only keep spikes that are more than STA_frame_num from the beginning
first_kept_frame = find(spike_times{1} > sta_frame_num, 1);
spike_times{1} = spike_times{1}(first_kept_frame:end);


%% Compute the STA
T=text_waitbar;

temp_sta = zeros(height, width, 3, sta_frame_num);

% put spike times into a vector for this RGC instead of a cell array
temp_spike_times = spike_times{1};
for fm = 1:length(temp_spike_times)
    
    T = text_waitbar(T,fm/length(temp_spike_times) - 0.01);   
    
    temp_mov = zeros(height,width,3,sta_frame_num);
    cntr = 0;
    for i= (temp_spike_times(fm) - sta_frame_num+1):temp_spike_times(fm)
        
        F = (mvi.getFrame(i-1).getBuffer - 0.5) * 2;
        
        F = reshape(F,3,width,height);
        cntr = cntr + 1;
        temp_mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end
    
    temp_sta = temp_sta + temp_mov;    
end

temp_sta = temp_sta ./ length(temp_spike_times);

%% inspect STA
for fm = 1:sta_frame_num
    image(norm_image(temp_sta(:,:,:,fm)))
    pause
end

%% calc RF, sig_stixels, and time_course from MATLAB STA
temp_rf = rf_from_sta(temp_sta);
params.marks_params.thresh = 4.5;
matlab_sig_stixels = significant_stixels(temp_sta, params.marks_params);
matlab_time_course = time_course_from_sta(temp_sta, sig_stixels);
figure(1)
image(norm_image(temp_rf))

figure(3)
matlab_tc = matlab_time_course(:,1) ./ norm(matlab_time_course(:,1));
plot(matlab_tc, 'k');
hold on

%% Check the acuracy of the MATLAB RF and TC against that produced in Vision
datarun = load_sta(datarun, 'load_sta', datarun.cell_ids(temp_cell_index));
vision_sta = datarun.stas.stas{temp_cell_index};
vision_sta = repmat(vision_sta, [1,1,3,1]);
vision_rf = rf_from_sta(vision_sta);
figure(2)
image(norm_image(vision_rf))
vision_sig_stixels = significant_stixels(vision_sta, params.marks_params);
vision_time_course = time_course_from_sta(vision_sta, sig_stixels);

%superimpose a red dashed line over the TC from matlab;
figure(3)
vision_tc = vision_time_course(16:30,1) ./ norm(vision_time_course(16:30,1));
plot(vision_tc, 'r--');




