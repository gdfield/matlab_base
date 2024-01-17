datapath = '/Users/gfield/Analysis/Flow/2018-09-25-0/data000/data000';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


%%
datapath_fl = '/Users/gfield/Analysis/Flow/2018-09-25-0/data001/data001';
datarun_fl = load_data(datapath_fl);
datarun_fl = load_neurons(datarun_fl);
datarun_fl = load_ei(datarun_fl, 'all', 'array_type', 519);

%%
datapath_fljt = '/Users/gfield/Analysis/Flow/2018-09-25-0/data002/data002';

datarun_fljt = load_data(datapath_fljt);
datarun_fljt = load_neurons(datarun_fljt);
datarun_fljt = load_ei(datarun_fljt, 'all', 'array_type', 519);


%% MAP to flow 
cell_types = {'ON brisk transient', 'ON brisk sustained', 'OFF brisk transient', 'OFF brisk sustained', 'OFF small transient', 'unclassified'};
[cell_map_to_flow, failed_cells] = map_ei(datarun, datarun_fl, 'master_cell_type', cell_types );

used_cell_indices = get_cell_indices(datarun, cell_types);
num_used_cells = length(used_cell_indices);

mapped_ids = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flow{used_cell_indices(rgc)})
        mapped_id(rgc) = 0;
    else
    mapped_ids(rgc) = cell_map_to_flow{used_cell_indices(rgc)};
    end
end

num_used_cells
num_mapped_cells = length(mapped_ids(mapped_ids >0))
fract_mapped_cells = num_mapped_cells/num_used_cells

mapped_ids_for_flow = mapped_ids(mapped_ids >0);

%% save out datarun and cell_id list
kill_fields = {'globals', 'piece', 'ei', 'channels'};
datarun_fl = rmfield(datarun_fl, kill_fields);

cd ~/Desktop
save datarun_fl datarun_fl
save mapped_ids_for_flow mapped_ids_for_flow


%% MAP to jittered flow
cell_types = {'ON brisk transient', 'ON brisk sustained', 'OFF brisk transient', 'OFF brisk sustained', 'OFF small transient', 'unclassified'};
[cell_map_to_flowjt, failed_cells] = map_ei(datarun, datarun_fljt, 'master_cell_type', cell_types );

used_cell_indices = get_cell_indices(datarun, cell_types);
num_used_cells = length(used_cell_indices);

mapped_ids = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flowjt{used_cell_indices(rgc)})
        mapped_id(rgc) = 0;
    else
    mapped_ids(rgc) = cell_map_to_flowjt{used_cell_indices(rgc)};
    end
end

num_used_cells
num_mapped_cells = length(mapped_ids(mapped_ids >0))
fract_mapped_cells = num_mapped_cells/num_used_cells

mapped_ids_for_flowjt = mapped_ids(mapped_ids >0);

%% save out datarun and cell_id list
kill_fields = {'globals', 'piece', 'ei', 'channels'};
datarun_fljt = rmfield(datarun_fljt, kill_fields);

cd ~/Desktop
save datarun_fljt datarun_fljt
save mapped_ids_for_flowjt mapped_ids_for_flowjt

%% Get the intersection of cells that mapped to both flow and jittered flow 

used_cell_indices = get_cell_indices(datarun, cell_types);
num_used_cells = length(used_cell_indices);
wn_cell_ids = datarun.cell_ids(used_cell_indices);

mapped_ids_flow = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flow{used_cell_indices(rgc)})
        mapped_ids_flow(rgc) = 0;
    else
    mapped_ids_flow(rgc) = cell_map_to_flow{used_cell_indices(rgc)};
    end
end

mapped_ids_flowjt = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flowjt{used_cell_indices(rgc)})
        mapped_ids_flowjt(rgc) = 0;
    else
    mapped_ids_flowjt(rgc) = cell_map_to_flowjt{used_cell_indices(rgc)};
    end
end

cell_maps = [wn_cell_ids', mapped_ids_flow, mapped_ids_flowjt];
save cell_maps cell_maps

%%
[cell_map_flow_to_flowjt, failed_cells] = map_ei(datarun_fl, datarun_fljt);

save cell_map_flow_to_flowjt cell_map_flow_to_flowjt

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New data set

datapath = '/Users/gfield/Analysis/Flow/2018-09-26-0/data007/data007';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);

%%
datapath_fl = '/Users/gfield/Analysis/Flow/2018-09-26-0/data005/data005';
datarun_fl = load_data(datapath_fl);
datarun_fl = load_neurons(datarun_fl);
datarun_fl = load_ei(datarun_fl, 'all', 'array_type', 519);

%%
datapath_fljt = '/Users/gfield/Analysis/Flow/2018-09-26-0/data006/data006';

datarun_fljt = load_data(datapath_fljt);
datarun_fljt = load_neurons(datarun_fljt);
datarun_fljt = load_ei(datarun_fljt, 'all', 'array_type', 519);

%% MAP to flow 
cell_types = {'ON brisk transient', 'ON brisk sustained', 'OFF brisk transient', 'OFF brisk sustained', 'OFF small transient', 'unclassified'};
[cell_map_to_flow, failed_cells] = map_ei(datarun, datarun_fl, 'master_cell_type', cell_types );
[cell_map_to_flowjt, failed_cells] = map_ei(datarun, datarun_fljt, 'master_cell_type', cell_types );


used_cell_indices = get_cell_indices(datarun, cell_types);
num_used_cells = length(used_cell_indices);
wn_cell_ids = datarun.cell_ids(used_cell_indices);

mapped_ids_flow = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flow{used_cell_indices(rgc)})
        mapped_ids_flow(rgc) = 0;
    else
    mapped_ids_flow(rgc) = cell_map_to_flow{used_cell_indices(rgc)};
    end
end

mapped_ids_flowjt = zeros(num_used_cells,1);
for rgc = 1:num_used_cells
    if isempty(cell_map_to_flowjt{used_cell_indices(rgc)})
        mapped_ids_flowjt(rgc) = 0;
    else
    mapped_ids_flowjt(rgc) = cell_map_to_flowjt{used_cell_indices(rgc)};
    end
end

cell_maps = [wn_cell_ids', mapped_ids_flow, mapped_ids_flowjt];
save cell_maps cell_maps

%% get the IDs of individual types:
ON_BS_IDs = get_cell_ids(datarun, 'ON brisk transient');
ON_BT_IDs = get_cell_ids(datarun, 'ON brisk sustained');
OFF_BS_IDs = get_cell_ids(datarun, 'OFF brisk transient');
OFF_BT_IDs = get_cell_ids(datarun, 'OFF brisk sustained');

save ON_BS_IDs ON_BS_IDs
save ON_BT_IDs ON_BT_IDs
save OFF_BS_IDs OFF_BS_IDs
save OFF_BT_IDs OFF_BT_IDs


%% save out datarun and cell_id list
kill_fields = {'globals', 'piece', 'ei', 'channels'};
datarun_fl = rmfield(datarun_fl, kill_fields);

cd ~/Desktop
save datarun_fl datarun_fl

%% save out datarun and cell_id list
kill_fields = {'globals', 'piece', 'ei', 'channels'};
datarun_fljt = rmfield(datarun_fljt, kill_fields);

cd ~/Desktop
save datarun_fljt datarun_fljt

%% save out datarun and cell_id list
kill_fields = {'globals', 'piece', 'ei', 'channels'};
datarun_wn = rmfield(datarun, kill_fields);

cd ~/Desktop
save datarun_wn datarun_wn


%% Get the intersection of cells that mapped to both flow and jittered flow 

%%
[cell_map_flow_to_flowjt, failed_cells] = map_ei(datarun_fl, datarun_fljt);

save cell_map_flow_to_flowjt cell_map_flow_to_flowjt





