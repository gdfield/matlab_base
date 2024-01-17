% load master data set
datapath = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data006/data006';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_sta(datarun, 'load_sta', 'all');

% get RFs and othe STA summaries.
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);
datarun = get_interspikeinterval(datarun, 'all');

%% load slave data set
datapath_s = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data005/data005';

datarun_s = load_data(datapath_s);
datarun_s = load_neurons(datarun_s);
datarun_s = load_params(datarun_s);
datarun_s = load_ei(datarun_s, 'all');
datarun_s = get_interspikeinterval(datarun_s, 'all');


%% plot and print as pdfs the rasters for chosen cell type

% choose cell type
temp_cell_type ='OFF brisk transient';

% map cells between master and slave datar un
[cell_list, failed_cells] = map_ei(datarun, datarun_s, 'master_cell_type', temp_cell_type);


% sort IDs and indices and get rid of cells that failed to map.
master_cell_indices = get_cell_indices(datarun, temp_cell_type);
mapped_list = cell_list(master_cell_indices);
% find empty fields
tp_ind = 0;
new_mapped_list = [];
new_master_indices = [];
for ct = 1:length(mapped_list)
    tmp_val = mapped_list{ct};
    if ~isempty(tmp_val)
        tp_ind = tp_ind+1;
        new_mapped_list(tp_ind) = mapped_list{ct};
        new_master_indices(tp_ind) = master_cell_indices(ct);
    end
end

% get indices of triggers that indicate beginning of each stimulus
epoch_trigger_inds = 1:6:length(datarun_s.triggers);
% get the actual triggers
epoch_triggers = datarun_s.triggers(epoch_trigger_inds);

% get indices in slave data
slave_indices = get_cell_indices(datarun_s, new_mapped_list);


figure(1); clf;
plot_axes = subplot_axes(1, [.07 .05 .93 .9],.1,.1,3,2);
orient(1,'landscape')
set(1,'defaultAxesFontSize',10)

% plot and print rasters
for ct = 1:length(slave_indices)

    % clear axes
    
    % plot the raster
    axes(plot_axes{1});cla
    get_raster(datarun_s.spikes{slave_indices(ct)}, epoch_triggers, 'stop', 9.9, 'foa', -1);
    title([temp_cell_type,' ',num2str(new_mapped_list(ct))])

    % plot ISI of slave
    axes(plot_axes{2});cla
    plot(datarun_s.interspikeinterval{slave_indices(ct)}.bins, datarun_s.interspikeinterval{slave_indices(ct)}.probabilities);
    set(gca,'xtick',[],'ytick',[]); 
    title('inter-spike intervals')
   
    % plot EI of slave
    axes(plot_axes{3});cla
    plot_ei(datarun_s, new_mapped_list(ct), 'foa', -1)
    set(gca,'xtick',[],'ytick',[]);  
    axis equal
    title('electrical image')
    
    % plot ISI of master
    axes(plot_axes{5});cla
    plot(datarun.interspikeinterval{new_master_indices(ct)}.bins, datarun.interspikeinterval{new_master_indices(ct)}.probabilities);
    %set(gca,'xtick',[],'ytick',[]);  
   
    % plot EI of master
    axes(plot_axes{6});cla
    plot_ei(datarun, datarun.cell_ids(new_master_indices(ct)), 'foa', -1)
    set(gca,'xtick',[],'ytick',[]);  
    axis equal

    axes(plot_axes{4});cla   
    plot_rf_summaries(datarun, temp_cell_type, 'plot_fits', true, 'fit_color', [0.5 0.5 0.5])
    plot_rf_summaries(datarun, datarun.cell_ids(new_master_indices(ct)), 'plot_fits', true, 'fit_color', [1 0 0], 'fit_width', 2, 'clear', false)
    set(gca,'xtick',[],'ytick',[]); 
    ylabel('master')
    
    file_name = ['~/Desktop/raster-plots/',num2str(new_mapped_list(ct)),'.pdf'];
    print(1, file_name, '-dpdf')
           
end
