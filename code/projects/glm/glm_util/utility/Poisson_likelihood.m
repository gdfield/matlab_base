function LL = Poisson_likelihood(firing_rate, recorded_raster, bindur)
% NJB 2016-06-03
% Takes in a firing rate and predicts Poisson spiking 
% The probability of a spike in each time bin is
% exp(-bin_duration*firing_rate)
% firing rate is in Hz (spikes per second)
% bin_duration is time (seconds) 
% bins_per_frame is the number of bins you should use for each entry of the
% firing rate 
% monitor refresh is in Hz (used to calculate the bin duration)
% trials is the number of trials you want to predict

if length(firing_rate) ==1
    firing_rate = firing_rate*ones(size(recorded_raster(1,:))); 
end
for i_trial = 1:size(recorded_raster,1)
    spikes = logical(recorded_raster(i_trial,:));
    nospikes = ~logical(recorded_raster(i_trial,:));
    %LL_trial(i_trial) = sum(log(bindur*firing_rate(spikes))) + sum(log(1-firing_rate(nospikes)*bindur));
    LL_trial(i_trial) = sum(log(bindur*firing_rate(spikes))) - sum(firing_rate(:)*bindur);
end
LL = sum(LL_trial);

end

