% extract the OSIs for all the cells in a chosen dataset


% define the data set to inspect
folderName = '2012-10-31-1';
dset = 5; % this needs to be updated to get that dataset in os_datasets.m


% get the paths to the data
data_list = os_datasets;
path_info = data_list(dset);

% load stuff grating data
datarun_g = load_data(path_info.grating_datapath);
datarun_g = load_neurons(datarun_g);
datarun_g = load_params(datarun_g);
datarun_g = load_ei(datarun_g, 'all');
datarun_g.names.stimulus_path = path_info.stimulus_path;
datarun_g = load_stim(datarun_g, 'user_defined_trigger_interval', path_info.trigger_interval);


[OSIs, DSIs, corr_vals] = find_OS_RGCs(datarun_g, 'all', 'plot_tuning', false,...
                            'plot_tuning_os_only', false, 'verbose', false);

%% Plot some example cell tuning curves

% extract spatial and temporal periods from datarun
spatial_periods = datarun_g.stimulus.params.SPATIAL_PERIOD;
temp_periods = datarun_g.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(datarun_g.stimulus.params.DIRECTION);

new_green = [78,178,101] ./ 256;
new_orange = [238,128,38] ./256;
new_blue = [97,149,207]./256;
new_gray = [0.5 0.5 0.5];

example_v_os_cell = 1549;
v_os_index = get_cell_indices(datarun_g, example_v_os_cell);
temp_spike_times = datarun_g.spikes{v_os_index};
fig_counter = 6;
for spat_p = 1:1
    for temp_p = 1:num_tps
        temp_sp = spatial_periods(spat_p);
        temp_tp = temp_periods(temp_p);
        fig_title = ['cell id-', num2str(datarun_g.cell_ids(v_os_index)), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
        [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_g.stimulus,...
                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
        plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_g.stimulus,...
                            'fig_num', fig_counter, 'fig_title', fig_title, 'raster_color', new_orange)  
        drawnow
        exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
        fig_counter = fig_counter + 1;
    end
end

example_h_os_cell = 60;
h_os_index = get_cell_indices(datarun_g, example_h_os_cell);
temp_spike_times = datarun_g.spikes{h_os_index};
fig_counter = 2;
for spat_p = 1:1
    for temp_p = 1:num_tps
        temp_sp = spatial_periods(spat_p);
        temp_tp = temp_periods(temp_p);
        fig_title = ['cell id-', num2str(datarun_g.cell_ids(h_os_index)), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
        [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_g.stimulus,...
                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
        plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_g.stimulus,...
                                'fig_num', fig_counter, 'fig_title', fig_title, 'raster_color', new_green)  
        drawnow
        exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
        fig_counter = fig_counter + 1;
    end
end

example_off_cell = 24;
off_cell_index = get_cell_indices(datarun_g, example_off_cell);
temp_spike_times = datarun_g.spikes{off_cell_index};
fig_counter = 10;
for spat_p = 1:1
    for temp_p = 1:num_tps
        temp_sp = spatial_periods(spat_p);
        temp_tp = temp_periods(temp_p);
        fig_title = ['cell id-', num2str(datarun_g.cell_ids(off_cell_index)), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
        [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_g.stimulus,...
                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
        plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_g.stimulus,...
                                'fig_num', fig_counter, 'fig_title', fig_title, 'raster_color', new_blue)  
        drawnow
        exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
        fig_counter = fig_counter + 1;
    end
end


%%
% calculat the mean DSI and mean response corr to remove DS cells and
% variabile response cells from OSI scatter plot.

num_rgcs = length(datarun_g.cell_ids);

% initialize some vectors to fill in the loop below
mean_DSI = zeros(num_rgcs,1);
mean_OSI = zeros(num_rgcs,1);
ave_corr = zeros(num_rgcs,1);

for rgc = 1:num_rgcs
    % get the mean DSI for each cell
    temp_DSI = DSIs(rgc, :,:);
    mean_DSI(rgc) = mean(temp_DSI(:));

    % get mean OSI for each cell
    temp_OSI = OSIs(rgc, :, :);
    mean_OSI(rgc) = mean(temp_OSI(:));

    temp_corr_strength = corr_vals(rgc, 1,:);
    ave_corr(rgc) = mean(temp_corr_strength(:));
end

% select thresholds for DSI on average response correlation
non_DSI_indices = find(mean_DSI < 0.3);
high_corr_indices = find(ave_corr > 0.8); 

OSIs_to_plot_indices = intersect(non_DSI_indices, high_corr_indices);

keep_flag = zeros(num_rgcs, 1);
keep_flag(OSIs_to_plot_indices) = 1;


figure(1); clf;
spat_per = [1 2];
temp_per = [2 1];
OSIs_to_plot = zeros(size(OSIs));
OSIs_to_plot(OSIs_to_plot_indices,:,:) = OSIs(OSIs_to_plot_indices, : ,:);
scatter(OSIs_to_plot(:, spat_per(1), temp_per(1)), OSIs_to_plot(:, spat_per(1), temp_per(2)), 'o', 'MarkerEdgeColor', new_gray, 'MarkerFaceColor', [1 1 1],'SizeData', 15)
hold on
scatter(OSIs_to_plot(h_os_index,2,1), OSIs_to_plot(h_os_index, 2,1), 'o', 'MarkerEdgeColor', new_green, 'MarkerFaceColor', new_green, 'SizeData', 40)
scatter(OSIs_to_plot(v_os_index,2,1), OSIs_to_plot(v_os_index, 2,1), 'o', 'MarkerEdgeColor', new_orange, 'MarkerFaceColor', new_orange, 'SizeData', 40)
scatter(OSIs_to_plot(off_cell_index,2,1), OSIs_to_plot(off_cell_index, 2,1), 'o', 'MarkerEdgeColor', new_blue, 'MarkerFaceColor', new_blue, 'SizeData', 40)
plot([0.3 0.3], [0 0.3], 'k--')
plot([0 0.3], [0.3 0.3], 'k--')

axis square
axis([0 0.6 0 0.6])
xticks([0 0.2 0.4 0.6])
yticks([0 0.2 0.4 0.6])
title('OSI')
xlabel('grating 1')
ylabel('grating 2')
exportgraphics(gca, 'osi_scatter.pdf', 'ContentType', 'vector')


%%

OSI_indices_g1 = find(OSIs(:, 1, 1) > 0.3);
OSI_indices_g2 = find(OSIs(:, 1, 2) > 0.3);

OSI_indices = intersect(OSI_indices_g1, OSI_indices_g2);

num_OSI_cells = length(OSI_indices);


for rgc = 1:num_OSI_cells
    for spat_p = 1:1
        for temp_p = 1:1
            temp_sp = spatial_periods(spat_p);
            temp_tp = temp_periods(temp_p);
            fig_title = ['cell id-', num2str(datarun_g.cell_ids(off_cell_index)), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
            [temp_tuning, temp_spike_nums, temp_rates] = get_direction_tuning(temp_spike_times, datarun_g.stimulus,...
                                                    'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
            plot_direction_tuning(temp_tuning, temp_spike_nums, datarun_g.stimulus,...
                                    'fig_num', fig_counter, 'fig_title', fig_title, 'raster_color', new_blue)  
            drawnow
            %exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
            fig_counter = fig_counter + 1;
        end
    end
end


