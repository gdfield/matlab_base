
% fft analysis of time courses

intervals = [4 2 1 1 1 1]; % frame refresh of the white noise stimulus
fs = 60./intervals;
samplerate = fs(j);

temp_tc = all_normd_tcs.(types_interest{type}).(all_ndf_names{j})(:,rgc); % grab the normalized time course (by the norm, not the peak)
[powerspec_xvalues, mean_powerspec, fft_signal] = PowerSpectrumFinder(temp_tc',samplerate) ;
freq_ranges.(all_ndf_names{j}) = powerspec_xvalues;
norm_ps = mean_powerspec/max(mean_powerspec); % normalize these values to the peak

% grab peak frequency
[maxs,inds] = max(norm_ps);
pk_freq = freq_ranges.(all_ndf_names{j})(inds);
