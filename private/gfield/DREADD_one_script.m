datapath = '/Users/gfield/Analysis/2020-11-04-0/data002/data002';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);


num_RGCs = length(datarun.cell_ids);

trig_indices = 1:6:length(datarun.triggers);
rep_triggers = datarun.triggers(trig_indices);

cd ~/Analysis/2020-11-04-0/

for rgc = 1:num_RGCs
    figure(1); clf;
    tmp_spike_times = datarun.spikes{rgc};
    temp_raster = get_raster(tmp_spike_times, rep_triggers, 'stop', 9.6, 'foa',1);
    drawnow
    tmp_name = ['Cell ', num2str(datarun.cell_ids(rgc))];
    title(tmp_name)
    orient(1,'landscape')
    print(1, tmp_name, '-dpdf')
end

