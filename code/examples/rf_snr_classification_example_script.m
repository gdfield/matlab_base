%% Spitting out deciles to Vision classification.txt
% datarun = load_data('/Volumes/Analysis/2016-06-13-1/data004/data004');
% tic; datarun = load_sta(datarun, 'save_rf', true); toc
% rf_snr_classification(datarun);


%% Getting the highest SNR cells from a class
piece = '2016-04-21-1';
data = '/data000/data000';
cellspec = [481 738 856 871 887 931 963 1023 1037 1083 1111 1113 1172 1201 1216 1231 1246 1261 1381 1396 1457 1472 1561 1651 1697 4246 4367 5088 ]; %OFF parasol

robust_std_method = 6;
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

datarun = load_data(fullfile(piece, data), loadopts);
datarun = load_sta(datarun, 'load_sta', []);
datarun = get_sta_summaries(datarun, 'all', 'robust_std_method', robust_std_method);
datarun = calc_rf_snrs(datarun);
snr =datarun.stas.medsigs(get_cell_indices(datarun, cellspec));
mean(snr)
figure; hist(snr);

% [sorted, inds] = sortlownans(datarun.stas.medsigs);
% snr_sorted_cells = datarun.cell_ids(inds(~isnan(sorted)));
% clear sorted inds
% 
% highest_snr_cells = snr_sorted_cells(end:-1:end-10)'