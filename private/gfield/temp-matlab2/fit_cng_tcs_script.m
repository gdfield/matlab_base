cd ~/Desktop
load datarun_5M

monitor_refresh = 60.35;

datarun = datarun_5M;
clear datarun_5M


for dset = 1:length(datarun)

    % extract indices of RGCs with space-time separable, high SNR STAs
    temp_rgc_indices = datarun{dset}.stas.keeper_indices;
    
    time_to_zero = zeros(length(temp_rgc_indices), 1);
    time_to_peak = zeros(length(temp_rgc_indices), 1);
    
    num_frames = length(datarun{dset}.stas.keeper_tcs(1,:));
    refresh_rate = datarun{dset}.stimulus.interval / monitor_refresh;
    time_bf_spike = -1*refresh_rate * (num_frames-1):refresh_rate:0;
    
    % loop over these rgcs to fit
    for rgc = 1:length(temp_rgc_indices)
        
        temp_tc = datarun{dset}.stas.keeper_tcs(rgc,:);

        %fit RGCs
        input_params.fit_n_one_filters = false;
        input_params.fit_n_two_filters = false;
        input_params.initial_n_one_filters = 9;
        input_params.initial_n_two_filters = 9;
        [check_fit_tc, fit_params] = fit_time_course(temp_tc', input_params);

        % pack fit parameters into structure and store in datarun
        tc_fit_params.scale_one = fit_params(1); %scale_one;
        tc_fit_params.scale_two = fit_params(2); %scale_two;
        tc_fit_params.tau_one = fit_params(3); %tau_one;
        tc_fit_params.tau_two = fit_params(4); %tau_two;
        tc_fit_params.n_one_filters = fit_params(5); %n_one_filters;
        tc_fit_params.n_two_filters = fit_params(6); %n_two_filters;

        datarun{dset}.matlab.sta_fits{temp_rgc_indices(rgc)} = tc_fit_params;

        extract_params.ms_per_frame = refresh_rate; 
        extract_params.time_range = 3:0.01:10;
        [temp_time_to_zero, temp_time_to_peak] = get_time_to_zero(datarun{dset}, datarun{dset}.cell_ids(temp_rgc_indices(rgc)), extract_params);

        % add zero crossing marker to indicate fit
        % plot as sanity check on fit
        figure(1); clf;
        plot(time_bf_spike,temp_tc, 'b')
        hold on
        plot(time_bf_spike,check_fit_tc, 'r')
        plot(time_bf_spike, zeros(length(time_bf_spike), 1), 'k')
        [dat_min, dat_max] = bounds([temp_tc, check_fit_tc]);
        plot([-1*temp_time_to_zero, -1*temp_time_to_zero], [dat_min, dat_max], 'k')
        legend('data', 'fit', 'baseline', 'zero cross', 'Location', 'NorthWest')
        title(num2str(datarun{dset}.cell_ids(temp_rgc_indices(rgc))))
        xlabel('s'); ylabel('amplitude'); axis square; axis tight        
        drawnow; pause(0.1)
        hold off

        
        time_to_zero(rgc) = temp_time_to_zero;
        time_to_peak(rgc) = temp_time_to_peak;
    end
    
    datarun{dset}.keeper_ttz = time_to_zero;
    datarun{dset}.keep_ttp = time_to_peak;
end

datarun_5M = datarun;
save datarun_5M datarun_5M

%%

all_dset_ttz = [];
for dset = 1:length(datarun)
    all_dset_ttz = [all_dset_ttz; datarun{dset}.keeper_ttz];
    all_dset_ttp = [all_dset_ttz; datarun{dset}.keep_ttp];
end

figure(5)
[h_min, h_max] = bounds(all_dset_ttz);
be_hist = linspace(h_min, h_max, 100);
[h_counts, ~] = histcounts(all_dset_ttz, be_hist);
histogram('BinCounts', h_counts, 'BinEdges', be_hist)


figure(6)
[h_min, h_max] = bounds(all_dset_ttp);
be_peaks = linspace(h_min, h_max, 100);
[h_peaks, ~] = histcounts(all_dset_ttp, be_peaks);
histogram('BinCounts', h_peaks, 'BinEdges', be_peaks)






