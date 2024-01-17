function [mean_tc, keeper_inds] = robust_tc_mean(tc_matrix, cor_threshold)

temp_mean = mean(tc_matrix);
tc_prjs = tc_matrix * temp_mean';

keeper_inds = find(tc_prjs > cor_threshold);

mean_tc = mean(tc_matrix(keeper_inds,:));


