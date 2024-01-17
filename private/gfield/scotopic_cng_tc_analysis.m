%% setup data paths for each time point
% set serverpath

serverpath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/';

%------%

% Month 1 Data
path_1m{1} = [serverpath, '2019-07-18-0/data005/data005'];
%path_1m{2} = [serverpath, '2019-07-19-0/data005/data005'];
%path_1m{3} = [serverpath, '2020-04-21-0/data006/data006'];

% data sets under scotopic condition
path_1m_child{1} = [serverpath, '2019-07-18-0/data003/data003'];
%path_1m_child{2} = [serverpath, '2019-07-19-0/data003/data003'];
%path_1m_child{3} = [serverpath, '2020-04-21-0/data003/data003'];


% loop over dataset to load and get_sta_summaries
opt_set = struct('load_neurons', true, 'load_sta', false, 'load_params', true);
for dset = 1:length(path_1m)
    dat_1M{dset} = load_data(path_1m{dset}, opt_set);
    dat_1M{dset} = load_ei(dat_1M{dset}, 'all', 'array_type', 512);
    
    dat_1M_child{dset} = load_data(path_1m_child{dset}, opt_set);
    dat_1M_child{dset} = load_ei(dat_1M_child{dset}, 'all', 'array_type', 512);
end


for dset = length(path_1m)
    [temp_map, ~] = map_ei(dat_1M{dset}, dat_1M_child{dset});
    map_to_child{1} = temp_map;
end

%------%



%% Analysis begins here!!!!
% Specify the time point to analyze
datarun_set = dat_1M;
child_set = path_1m_child;

% choose cell types to analyze
cell_types = {'ON brisk transient'};

% set threshold parameter for getting significant stixels in get_sta_summaries
marks_params.thresh = 3.5;

% loop through cell types and dataruns to accumulate time courses
clear time_course_struct
for ct = 1:length(cell_types)
time_courses = [];
    
    for dset = 1:length(datarun_set)
        
        % get IDs in the child data set for specified cell type.
        
        temp_cell_indices = get_cell_indices(datarun_set{dset}, cell_types);
        tmp_array = map_to_child{dset};
        for cc = 1:length(temp_cell_indices)
           
            child_dset_IDs = tmp_array{temp_cell_indices(cc)}
        end
        %remove empty fields
        child_dset
        
        datarun_set{dset} = get_sta_summaries(datarun_set{dset}, child_dset_IDs, 'marks_params', marks_params);
        temp_tcs = get_time_courses_matrix(datarun_set{dset}, cell_types{ct});
        time_course_struct{ct}.time_courses{dset}.time_c = temp_tcs;
        %time_courses = [time_courses; temp_tcs];
        %x = size(time_courses)
    end    

    time_course_struct{ct}.cell_type = cell_types{ct};
%     time_course_struct{ct}.time_courses = time_courses;
    
end


fm_interval = dat_1M{1}.stimulus.interval;
num_frames = size(dat_1M{1}.stas.stas{1},4);
refresh = 1000/60.35;
frame_time = fm_interval * refresh;
time_points = -1*(num_frames-1)*frame_time:frame_time:0;


%% plot time courses
mouse_age = input('At what time point was this data collected? ', 's');
for ct = 1:length(cell_types)
    figure(ct); clf
    plot(time_points, zeros(length(time_points),1), 'k') % plot zero line
    hold on
    for tm = 1:length(time_course_struct{ct}.time_courses)
        plot(time_points, time_course_struct{ct}.time_courses{tm}.time_c') % plot data
        hold on
    end
    axis square
    axis tight
    title([mouse_age, ': ', cell_types{ct}])
    xlabel('time (ms)')
    ylabel('amplitude (a.u.)')
end

%% STA SNR analysis

clear spike_counts rf_snr
for ct = 1:length(cell_types)
    
    for dset = 1:length(datarun_set)
        
        temp_cell_ids = get_cell_indices(datarun_set{dset}, cell_types{ct});
        temp_num_rgcs = length(temp_cell_ids);
        temp_num_rates = zeros(temp_num_rgcs, 1);
        snrs = zeros(temp_num_rgcs, 1);
        
        for rgc = 1:temp_num_rgcs
            
            temp_num_rates(rgc) = length(datarun_set{dset}.spikes{temp_cell_ids(rgc)}) ./ datarun_set{dset}.duration; %datarun_set{dset}.spikes{temp_cell_ids(rgc)}
            temp_marks = datarun_set{dset}.stas.marks{temp_cell_ids(rgc)};
            %temp_marks is empty - rest of code is NaN or 0
            temp_rf = datarun_set{dset}.stas.rfs{temp_cell_ids(rgc)};
            temp_rf = temp_rf ./ norm(temp_rf(:));
            temp_vals = temp_rf(temp_marks);
            temp_sig = mean(temp_vals);
            temp_noise = robust_std(temp_rf(:));
            temp_snr = temp_sig ./ temp_noise;
        
            snrs(rgc) = temp_snr;
            
        end
        
        spike_counts(ct,dset).mean = mean(temp_num_rates);
        spike_counts(ct,dset).sd = std(temp_num_rates);
        spike_counts(ct,dset).raw_vals = temp_num_rates;
        
        rf_snr(ct,dset).mean = mean(snrs);
        rf_snr(ct,dset).sd = std(snrs);
        rf_snr(ct,dset).snrs = snrs;        
        
    end
end
        

figure(2); clf;
errorbar(rf_snr(1,1).mean, spike_counts(1, 1).mean, spike_counts(1, 1).sd, 'ko')
hold on
errorbar(rf_snr(1,2).mean, spike_counts(1, 2).mean, spike_counts(1, 2).sd,'k^')
errorbar(rf_snr(1,3).mean, spike_counts(1, 2).mean, spike_counts(1, 3).sd,'k^')
errorbar(rf_snr(2,1).mean, spike_counts(2, 1).mean, spike_counts(2, 1).sd,'ro')
errorbar(rf_snr(2,2).mean, spike_counts(2, 2).mean, spike_counts(2, 2).sd,'r^')
errorbar(rf_snr(2,3).mean, spike_counts(2, 3).mean, spike_counts(2, 2).sd,'r*')


plot(rf_snr(1,1).snrs, spike_counts(1,1).raw_vals, 'k.')
plot(rf_snr(1,2).snrs, spike_counts(1,2).raw_vals, 'k.')
plot(rf_snr(1,3).snrs, spike_counts(1,2).raw_vals, 'k.')
plot(rf_snr(2,1).snrs, spike_counts(2,1).raw_vals, 'r.')
plot(rf_snr(2,2).snrs, spike_counts(2,2).raw_vals, 'r.')
plot(rf_snr(2,3).snrs, spike_counts(2,2).raw_vals, 'r.')

xlabel('SNR of STA');
ylabel('spike rate');
text(16, 42, 'red = OFF BT; black = ON BT')
hold off
    

%% compute and fit the static nonlinearities

% Get SNLs and fits and save to disk.

for dset = 1:length(datarun_set)
    datarun_set{dset} = load_java_movie(datarun_set{dset}, '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-2-0.48-11111-52x40-60.35.xml');
    datarun_set{dset} = get_snls(datarun_set{dset}, cell_types, 'verbose', true, 'new', true, 'fit', 'cum_norm');
end

% save SNL info
cd ~/Desktop
RGC_SNLs = datarun.stas.snls;
save RGC_SNLs RGC_SNLs

temp_indices = get_cell_indices(datarun, {'ON type1', 'ON type2', 'ON type3','OFF type1', 'OFF type2', 'OFF type3'});
num_rgcs = length(temp_indices);









