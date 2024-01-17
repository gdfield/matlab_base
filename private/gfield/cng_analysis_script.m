%% analyze spontaneous firing

%% load data from WN run
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-10-0/data000/data000';

datarun_wn = load_data(wn_datapath);
datarun_wn = load_neurons(datarun_wn);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_ei(datarun_wn, 'all');


%% load data from spontaneous activity

spont_dark_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-10-0/data000/data000';

datarun = load_data(spont_dark_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

%% map from WN runs to spontaneous activity

%define cell type
cell_type = 'all';
corr_threshold = 0.95;

% map from wn to spontaneous run using EIs
[map_list, failed_maps] = map_ei(datarun_wn, datarun, 'master_cell_type', cell_type, 'corr_threshold', corr_threshold);

% extract cell_ids from the CELL ARRAY map_list
tmp_cntr = 0;
clear slave_ids;
for rgc = 1:length(map_list)
    if ~isempty(map_list{rgc})
        tmp_cntr = tmp_cntr + 1;
        slave_ids(tmp_cntr)= map_list{rgc};      
    end
end


%% Compute the mean firing rate
temp_cell_indices = get_cell_indices(datarun, slave_ids);

% calculate mean spike rate
spike_rates = zeros(length(temp_cell_indices), 1);
for rgc = 1:length(temp_cell_indices)
    temp_num_spikes = length(datarun.spikes{temp_cell_indices(rgc)});
    temp_rate = temp_num_spikes ./ datarun.duration;
    spike_rates(rgc) = temp_rate;
end

mean_spike_rate = mean(spike_rates); % Hz


%% Compute the power spectrum

%cell_type = 'all';
cell_type = [1, 16, 31, 46];
temp_cell_indices = get_cell_indices(datarun, cell_type);

num_samples = datarun.duration * datarun.sampling_rate;


% calculate mean spike rate
ps_of_spike_trains = zeros(length(temp_cell_indices),num_samples/2);
for rgc = 1:length(temp_cell_indices)

    % get spike train and convert to units of samples
    temp_spike_train = floor(datarun.spikes{temp_cell_indices(rgc)} * datarun.sampling_rate);

    % translate spike times into a string of zeros and ones at the sampled bin size
    binned_spike_train = zeros(num_samples, 1);
    binned_spike_train(temp_spike_train) = 1;

    % convert the duration and sampling rate into a frequency range
    temp_freq = [0:num_samples/2-1] / datarun.duration;

    % compute the power spectrum (ps) from the fft
    ps_spike_train = abs(fft(binned_spike_train)) ./ (num_samples./2);
    ps_spike_train = ps_spike_train(1:num_samples./2).^2;
    ps_spike_train = 10*log10(ps_spike_train); % convert to dB

    ps_of_spike_trains(rgc, :) = ps_spike_train;
end


% inspect
plot(temp_freq, ps_of_spike_trains(1,:))
axis([0 35 -200 0])
label('power spectrum')
xlabel('frequency (Hz)')
ylabel('power db/Hz')

