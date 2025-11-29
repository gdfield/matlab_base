%% Plot the mosaic of the h-OS cells

datapath = '/Volumes/gdf/rat-data/2012-10-31-0/YASS/data000/data000';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'OFF h os');

marks_params.thresh = 3.5;
datarun = get_sta_summaries(datarun, 'OFF h os', 'marks_params', marks_params);

% set saving path
cd /Users/gfield/Desktop/os_analysis/2012-10-31-0/rf_plots

% plot h os cells w/o labels
plot_rf_summaries(datarun, 'OFF h os', 'plot_fit', true)
title('horizontal OS simple cells')
exportgraphics(gcf, 'hos_mosaic.pdf', 'ContentType', 'vector')

plot_rf_summaries(datarun, 'OFF h os', 'plot_fit', true, 'label', true)
title('horizontal OS simple cells')
exportgraphics(gcf, 'hos_mosaic_labels.pdf', 'ContentType', 'vector')


%% Plot the RFs of all of the definite OS cells

% load white noise run
datapath = '/Volumes/gdf/rat-data/2012-10-31-0/YASS/data000/data000';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_ei(datarun, 'all', 'array_type', 512);

% set some parameters and get sta summaries
marks_params.thresh = 3.5;
filt_rad = 0.75;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);
datarun = set_polarities(datarun, 'guess', true);

% load gratings run
datapath_dg = '/Volumes/gdf/rat-data/2012-10-31-0/YASS/data002/data002';
datarun_dg = load_data(datapath_dg);
datarun_dg = load_neurons(datarun_dg);
datarun_dg = load_params(datarun_dg);
datarun_dg = load_ei(datarun_dg, 'all', 'array_type', 512);

% load the os list and maybe-list info
load /Users/gfield/Desktop/os_analysis/2012-10-31-0_os.mat

% may the os cell IDs from gratings to white noise run
mapped_cell_list = map_ei(datarun_dg, datarun, 'master_cell_type', os_cell_list, 'corr_threshold', 0.9);
os_cell_indices = get_cell_indices(datarun_dg, os_cell_list);

% ensure the polarity of the RF plots are correct
temp_tcs = get_time_courses_matrix(datarun, 'OFF h os', 'norm_flag', true);
tc_ext_val = ext(mean(temp_tcs));
if tc_ext_val > 0
    polarity_enforcer = 1;
elseif tc_ext_val < 0
    polarity_enforcer = -1;
end


% plot the RFs
cd /Users/gfield/Desktop/os_analysis/2012-10-31-0/rf_plots
for rgc = 1:length(os_cell_list)

    wn_cell_id = mapped_cell_list{os_cell_indices(rgc)};
    if isempty(wn_cell_id)
        warning('white noise cell ID is empty')
        continue
    end
    wn_cell_index = get_cell_indices(datarun, wn_cell_id);

    temp_rf = datarun.stas.rfs{wn_cell_index};

    temp_tc = datarun.stas.time_courses{wn_cell_index};
    tc_ext_val = ext(temp_tc);
    if tc_ext_val > 0
        polarity_enforcer = 1;
    elseif tc_ext_val < 0
        polarity_enforcer = -1;
    end


    temp_params = struct('filt_type','gauss','radius',filt_rad);
    [filt_rf, ~] = rf_filtered(temp_rf, temp_params);
    norm_rf = norm_image(polarity_enforcer * filt_rf);
    imagesc(squeeze(norm_rf(:,:,1)))
    colormap(brewermap([],'RdBu'))
    caxis([0,1])
    axis equal
    axis tight
    plot_rf_summaries(datarun, wn_cell_id, 'foa', 1, 'plot_fits', true)
    temp_title = ['wn', num2str(wn_cell_id),'-dg', num2str(os_cell_list(rgc))];
    title(temp_title)
    exportgraphics(gca, [temp_title, '.pdf'], 'ContentType','vector')

end

%% Plot the Time courses of the h OS cells
% plot timecourses
% compute temporal axis

datapath = '/Volumes/gdf/2012-10-31-0/YASS/data000/data000';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_ei(datarun, 'all', 'array_type', 512);
marks_params.thresh = 3.5;

datarun = get_sta_summaries(datarun, 'OFF h os', 'marks_params', marks_params);

frame_refresh = datarun.stimulus.interval;
monitor_refresh = datarun.stimulus.monitor_refresh;
STA_length = size(datarun.stas.stas{1}, 4);
base_dwell_time = 1000/monitor_refresh; % units of ms
stim_dwell_time = base_dwell_time * frame_refresh;
begin_time = (STA_length-1) * stim_dwell_time;
time_bins = 0:stim_dwell_time:begin_time;
temp_tcs = get_time_courses_matrix(datarun, 'OFF h os', 'norm_flag', true);
plot(time_bins, temp_tcs', 'k')
xlabel('frames')
ylabel('a.u.')
title('OFF TCs')
axis tight
exportgraphics(gca, 'off_tcs.pdf', 'ContentType','vector')


temp_indices = get_cell_indices(datarun, 'OFF h os');
temp_marks = significant_stixels(datarun.stas.stas{temp_indices(1)}, 'select', 'max');
imagesc(datarun.stas.rfs{temp_indices(1)})

imagesc(full(temp_marks))

%% Plot the RFs of all of the maybe OS cells

cd /Users/gfield/Desktop/os_analysis/2012-10-31-0/maybe_rf_plots

mapped_cell_list = map_ei(datarun_dg, datarun, 'master_cell_type', os_maybe_list, 'corr_threshold', 0.9);
os_cell_indices = get_cell_indices(datarun_dg, os_maybe_list);

for rgc = 1:length(os_maybe_list)

    wn_cell_id = mapped_cell_list{os_cell_indices(rgc)};
    if isempty(wn_cell_id)
        warning('white noise cell ID is empty')
        continue
    end
    wn_cell_index = get_cell_indices(datarun, wn_cell_id);

    temp_rf = datarun.stas.rfs{wn_cell_index};

    temp_tc = datarun.stas.time_courses{wn_cell_index};
    tc_ext_val = ext(temp_tc);
    if tc_ext_val > 0
        polarity_enforcer = 1;
    elseif tc_ext_val < 0
        polarity_enforcer = -1;
    end

    temp_params = struct('filt_type','gauss','radius',filt_rad);
    [filt_rf, ~] = rf_filtered(temp_rf, temp_params);
    norm_rf = norm_image(filt_rf);
    imagesc(squeeze(norm_rf(:,:,1)))
    colormap(brewermap([],'RdBu'))
    caxis([0,1])
    axis equal
    axis tight
    plot_rf_summaries(datarun, wn_cell_id, 'foa', 1, 'plot_fits', true)

    temp_title = ['wn', num2str(wn_cell_id),'-dg', num2str(os_maybe_list(rgc))];
    title(temp_title)
    exportgraphics(gca, [temp_title, '.pdf'], 'ContentType','vector')

end

%% Plot an OFF mosaic

% choose cell type
temp_cell_type = 'OFF type2';
datarun = set_polarities(datarun, 'guess', true);
filt_rad = 0.75;

% get STA summaries
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, temp_cell_type, 'marks_params', marks_params);

% set path to save stuff
cd /Users/gfield/Desktop/os_analysis/2012-10-31-0/off_rf_plots

% plot the mosaic
plot_rf_summaries(datarun, temp_cell_type, 'plot_fit', true)
title('OFF brisk transient RGCs')
exportgraphics(gcf, 'off_mosaic.pdf', 'ContentType', 'vector')

plot_rf_summaries(datarun, temp_cell_type, 'plot_fit', true, 'label', true)
title('OFF brisk transient RGCs')
exportgraphics(gcf, 'off_mosaic_labels.pdf', 'ContentType', 'vector')

% ensure the polarity of the RF plots are correct
temp_tcs = get_time_courses_matrix(datarun, temp_cell_type, 'norm_flag', true);
tc_ext_val = ext(mean(temp_tcs));
if tc_ext_val > 0
    polarity_enforcer = 1;
elseif tc_ext_val < 0
    polarity_enforcer = -1;
end

% plot the RFs
temp_cell_indices = get_cell_indices(datarun, temp_cell_type);
for rgc = 1:length(temp_cell_indices)

    temp_rf = datarun.stas.rfs{temp_cell_indices(rgc)};
    temp_params = struct('filt_type','gauss','radius',filt_rad);
    [filt_rf, ~] = rf_filtered(temp_rf, temp_params);
    norm_rf = norm_image(polarity_enforcer * filt_rf);
    imagesc(squeeze(norm_rf(:,:,1)))
    colormap(brewermap([],'RdBu'))
    caxis([0,1])
    axis equal
    axis tight
    plot_rf_summaries(datarun, datarun.cell_ids(temp_cell_indices(rgc)), 'foa', 1, 'plot_fits', true)
    temp_title = num2str(datarun.cell_ids(temp_cell_indices(rgc)));
    title(temp_title)
    exportgraphics(gca, [temp_title, '.pdf'], 'ContentType','vector')

end

% plot timecourses
% compute temporal axis
frame_refresh = datarun.stimulus.interval;
monitor_refresh = datarun.stimulus.monitor_refresh;
STA_length = size(datarun.stas.stas{1}, 4);
base_dwell_time = 1000/monitor_refresh; % units of ms
stim_dwell_time = base_dwell_time * frame_refresh;
begin_time = (STA_length-1) * stim_dwell_time;
end_time = 0;
time_bins = 0:stim_dwell_time:begin_time
temp_tcs = get_time_courses_matrix(datarun, temp_cell_type, 'norm_flag', true);
plot(time_bins, temp_tcs', 'k')
xlabel('frames')
ylabel('a.u.')
title('OFF TCs')
axis tight
exportgraphics(gca, 'off_tcs.pdf', 'ContentType','vector')

%% Plot an ON mosaic

% choose cell type
temp_cell_type = 'ON type2';
datarun = set_polarities(datarun, 'guess', true);
filt_rad = 0.75;

% get STA summaries
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, temp_cell_type, 'marks_params', marks_params);

% set path to save stuff
cd /Users/gfield/Desktop/os_analysis/2012-10-31-0/on_rf_plots

% plot the mosaic
plot_rf_summaries(datarun, temp_cell_type, 'plot_fit', true)
title('OFF brisk transient RGCs')
exportgraphics(gcf, 'off_mosaic.pdf', 'ContentType', 'vector')

plot_rf_summaries(datarun, temp_cell_type, 'plot_fit', true, 'label', true)
title('OFF brisk transient RGCs')
exportgraphics(gcf, 'off_mosaic_labels.pdf', 'ContentType', 'vector')


% ensure the polarity of the RF plots are correct
temp_tcs = get_time_courses_matrix(datarun, temp_cell_type, 'norm_flag', true);
tc_ext_val = ext(mean(temp_tcs));
if tc_ext_val > 0
    polarity_enforcer = 1;
elseif tc_ext_val < 0
    polarity_enforcer = -1;
end

% plot the RFs
temp_cell_indices = get_cell_indices(datarun, temp_cell_type);
for rgc = 1:length(temp_cell_indices)

    temp_rf = datarun.stas.rfs{temp_cell_indices(rgc)};
    temp_params = struct('filt_type','gauss','radius',filt_rad);
    [filt_rf, ~] = rf_filtered(temp_rf, temp_params);
    norm_rf = norm_image(polarity_enforcer * filt_rf);
    imagesc(squeeze(norm_rf(:,:,1)))
    colormap(brewermap([],'RdBu'))
    caxis([0,1])
    axis equal
    axis tight
    plot_rf_summaries(datarun, datarun.cell_ids(temp_cell_indices(rgc)), 'foa', 1, 'plot_fits', true)
    temp_title = num2str(datarun.cell_ids(temp_cell_indices(rgc)));
    title(temp_title)
    exportgraphics(gca, [temp_title, '.pdf'], 'ContentType','vector')

end

% plot timecourses
% compute temporal axis
frame_refresh = datarun.stimulus.interval;
monitor_refresh = datarun.stimulus.monitor_refresh;
STA_length = size(datarun.stas.stas{1}, 4);
base_dwell_time = 1000/monitor_refresh; % units of ms
stim_dwell_time = base_dwell_time * frame_refresh;
begin_time = (STA_length-1) * stim_dwell_time;
end_time = 0;
time_bins = 0:stim_dwell_time:begin_time
temp_tcs = get_time_courses_matrix(datarun, temp_cell_type, 'norm_flag', true);
plot(time_bins, temp_tcs', 'k')
xlabel('frames')
ylabel('a.u.')
title('OFF TCs')
axis tight
exportgraphics(gca, 'on_tcs.pdf', 'ContentType','vector')


