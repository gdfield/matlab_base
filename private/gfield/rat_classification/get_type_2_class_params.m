function class_params = get_type_2_class_params(datarun_wn, datarun_dg, cell_ids_wn, cell_ids_dg)

epoch_duration = 8;
num_epochs = length(datarun_dg.epoch_triggers_256);
% check that cell_ids lengths are the same
if cell_ids_wn ~= cell_ids_dg
    error('cell_ids for the white noise and grating conditions are now equal in length')
end

cell_indices_wn = get_cell_indices(datarun_wn, cell_ids_wn);
cell_indices_dg = get_cell_indices(datarun_dg, cell_ids_dg);

phasic_indices = zeros(length(cell_indices_wn),1);
spike_cntr_256 = zeros(length(cell_indices_wn),1);
spike_cntr_512 = zeros(length(cell_indices_wn),1);
wn_spike_rate = zeros(length(cell_indices_wn),1);
tf_to_pk = zeros(length(cell_indices_wn),1);
for rgc = 1:length(cell_indices_wn)
    % phasic index
    tmp_tc = datarun_wn.stas.time_courses{cell_indices_wn(rgc)};
    phasic_indices(rgc) = abs(sum(tmp_tc) ./ sum(abs(tmp_tc)));
    tmp_tc = tmp_tc ./ ext(tmp_tc);
    tf_to_pk(rgc) = min(tmp_tc) ./ max(tmp_tc);
    
    % mean response to low TP grating
    tmp_raster = get_raster(datarun_dg.spikes{cell_indices_dg(rgc)}, datarun_dg.epoch_triggers_256, 'stop', epoch_duration, 'plot', false);
    tmp_cntr = 0;
    for epch = 1:length(tmp_raster)
        tmp_cntr = tmp_cntr + length(tmp_raster{epch});
    end
    spike_cntr_256(rgc) = tmp_cntr ./ (epoch_duration * num_epochs);
    
    % mean response to high TP grating
    tmp_raster = get_raster(datarun_dg.spikes{cell_indices_dg(rgc)}, datarun_dg.epoch_triggers_512, 'stop', epoch_duration, 'plot', false);
    tmp_cntr = 0;
    for epch = 1:length(tmp_raster)
        tmp_cntr = tmp_cntr + length(tmp_raster{epch});
    end
    spike_cntr_512(rgc) = tmp_cntr ./ (epoch_duration * num_epochs);
     
    tmp_num_spikes = length(datarun_wn.spikes{cell_indices_wn(rgc)}) ./ datarun_wn.duration;  
    wn_spike_rate(rgc) = tmp_num_spikes; 
    
end

class_params.phasic_indices = phasic_indices;
class_params.spike_rate_256 = spike_cntr_256;
class_params.spike_rate_512 = spike_cntr_512;
class_params.wn_spike_rate = wn_spike_rate;
class_params.tf_to_pk = tf_to_pk;