datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data006/data006';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun,'load_sta', 'all');
datarun = load_ei(datarun, 'all');


datapath_rps = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data005/data005';
datarun_rps = load_data(datapath_rps);
datarun_rps = load_neurons(datarun_rps);
datarun_rps = load_params(datarun_rps);
datarun_rps = load_ei(datarun_rps, 'all');

datapath_nat = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data008/data008';
datarun_nat = load_data(datapath_nat);
datarun_nat = load_neurons(datarun_nat);
datarun_nat = load_params(datarun_nat);
datarun_nat = load_ei(datarun_nat, 'all');


%% mapping tp WN

cell_type = 'OFF brisk transient';
[cell_list, failed_list] = map_ei(datarun, datarun_rps, 'master_cell_type', cell_type);
temp_indices = get_cell_indices(datarun, cell_type);
new_list = cell_list(temp_indices);

num_mapped_cells = length(temp_indices) - length(failed_list);
wn_rep_matches = zeros(num_mapped_cells, 2);
cntr = 0;
for rgc = 1:length(temp_indices)
    if ~isempty(cell_list{temp_indices(rgc)})
        cntr = cntr +1;
        wn_rep_matches(cntr, 1) = datarun.cell_ids(temp_indices(rgc));
        wn_rep_matches(cntr, 2) = cell_list{temp_indices(rgc)};
    end
end

OFFbt_wn_rep_matches = wn_rep_matches;
save OFFbt_wn_rep_matches OFFbt_wn_rep_matches

%% mapping to nat movie

cell_type = 'OFF brisk transient';
[cell_list, failed_list] = map_ei(datarun, datarun_nat, 'master_cell_type', cell_type);
temp_indices = get_cell_indices(datarun, cell_type);
new_list = cell_list(temp_indices);

num_mapped_cells = length(temp_indices) - length(failed_list);
nat_rep_matches = zeros(num_mapped_cells, 2);
cntr = 0;
for rgc = 1:length(temp_indices)
    if ~isempty(cell_list{temp_indices(rgc)})
        cntr = cntr +1;
        nat_rep_matches(cntr, 1) = datarun.cell_ids(temp_indices(rgc));
        nat_rep_matches(cntr, 2) = cell_list{temp_indices(rgc)};
    end
end

OFFbt_nat_rep_matches = nat_rep_matches;
save OFFbt_nat_rep_matches OFFbt_nat_rep_matches
