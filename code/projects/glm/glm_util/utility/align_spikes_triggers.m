function spikes_adj = align_spikes_triggers(spikes, triggers, frames_per_trigger)
% NB 2016-06-3
% inputs: 
% spikes: spike times in seconds
% triggers: datarun.triggers (the trigger times in seconds
% frames_per_trigger: pretty self explanatory I'd say

% This function aligns the spikes to the triggers. Basically it calculates
% where the trigger should've been in an ideal world, and then moves the
% according to the difference between that time and the actual trigger
% time.
    monitor_refresh = frames_per_trigger/median(diff(triggers)); % calibrate monitor_refresh from triggers and frames per trigger 
    spikes_adj = [];
    n_block=0;
    block_length = frames_per_trigger/monitor_refresh;
    for i=1:(length(triggers))
        actual_t_start=triggers(i); % When was the trigger?
        supposed_t_start=n_block*block_length; % When was the trigger supposed to be?
        idx1=spikes > actual_t_start; % Find the spikes in this chunk of data
        idx2=spikes < (actual_t_start+block_length);
        spikes_adj=[spikes_adj; spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start]; % move those spikes by the difference between the ideal and actual trigger times
        n_block=n_block+1; % move on to the next chunk of time!
    end
end