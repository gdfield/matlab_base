PSAM_cell_indices = get_cell_indices(datarun, mapped_type_1);

% parse the PSAM cell indices for the drug and wash conditions
drug_cell_ids = zeros(1,length(PSAM_cell_indices));
wash_cell_ids = zeros(1,length(PSAM_cell_indices));
for cc = 1:length(PSAM_cell_indices)
    drug_cell_ids(cc) = drug_list_map{PSAM_cell_indices(cc)};
    wash_cell_ids(cc) = wash_list_map{PSAM_cell_indices(cc)};
end
num_PSAM_cells = length(PSAM_cell_indices);

% get means and SEs for pre, drug and wash
pre_mean_tc = average_time_course(datarun, mapped_type_1, 'norm_one', false, 'peak_norm', false);
pre_se_tc = standard_error_time_course(datarun, mapped_type_1, 'norm_one', false, 'peak_norm', false);
drug_mean_tc = average_time_course(drug_datarun, drug_cell_ids, 'norm_one', false, 'peak_norm', false);
drug_se_tc = standard_error_time_course(drug_datarun, drug_cell_ids, 'norm_one', false, 'peak_norm', false);
wash_mean_tc = average_time_course(wash_datarun, wash_cell_ids, 'norm_one', false, 'peak_norm', false);
wash_se_tc = standard_error_time_course(wash_datarun, wash_cell_ids, 'norm_one', false, 'peak_norm', false);

shadedErrorBar([],pre_mean_tc, pre_se_tc*1.96, 'k');
hold on
shadedErrorBar([],drug_mean_tc, drug_se_tc*1.96, 'r');
shadedErrorBar([],wash_mean_tc, wash_se_tc*1.96, 'b');
hold off
%axis([1 30 -0.7 0.6])
axis square

print(1, '~/Desktop/ave-tcs-PSAM-sensitive.pdf', '-dpdf')

%%
PSAM_cell_indices = get_cell_indices(datarun, mapped_type_2);

% parse the PSAM cell indices for the drug and wash conditions
drug_cell_ids = zeros(1,length(PSAM_cell_indices));
wash_cell_ids = zeros(1,length(PSAM_cell_indices));
for cc = 1:length(PSAM_cell_indices)
    drug_cell_ids(cc) = drug_list_map{PSAM_cell_indices(cc)};
    wash_cell_ids(cc) = wash_list_map{PSAM_cell_indices(cc)};
end
num_PSAM_cells = length(PSAM_cell_indices);

% get means and SEs for pre, drug and wash
pre_mean_tc = average_time_course(datarun, mapped_type_2, 'norm_one', true, 'peak_norm', false);
pre_se_tc = standard_error_time_course(datarun, mapped_type_2, 'norm_one', true, 'peak_norm', false);
drug_mean_tc = average_time_course(drug_datarun, drug_cell_ids, 'norm_one', true, 'peak_norm', false);
drug_se_tc = standard_error_time_course(drug_datarun, drug_cell_ids, 'norm_one', true, 'peak_norm', false);
wash_mean_tc = average_time_course(wash_datarun, wash_cell_ids, 'norm_one', true, 'peak_norm', false);
wash_se_tc = standard_error_time_course(wash_datarun, wash_cell_ids, 'norm_one', true, 'peak_norm', false);

shadedErrorBar([],pre_mean_tc, pre_se_tc*1.96, 'k');
hold on
shadedErrorBar([],drug_mean_tc, drug_se_tc*1.96, 'r');
shadedErrorBar([],wash_mean_tc, wash_se_tc*1.96, 'b');
hold off
%axis([1 30 -0.7 0.6])
axis square

print(1, '~/Desktop/ave-tcs-PSAM-sensitive.pdf', '-dpdf')
