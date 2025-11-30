function subsampled_times = subsample_triggers(trigger_times, min_spacing)
% SUBSAMPLE_TRIGGERS Subsample trigger times with minimum spacing
%
% Inputs:
%   trigger_times - vector of trigger times (frame indices)
%   min_spacing - minimum number of frames between triggers (e.g., 50)
%
% Outputs:
%   subsampled_times - subsampled vector of trigger times

    % Sort trigger times (in case they're not already sorted)
    trigger_times = sort(trigger_times(:));
    
    if isempty(trigger_times)
        subsampled_times = [];
        return;
    end
    
    % Initialize with first trigger
    subsampled_times = trigger_times(1);
    
    % Greedily select triggers that are at least min_spacing away
    for i = 2:length(trigger_times)
        if trigger_times(i) - subsampled_times(end) >= min_spacing
            subsampled_times = [subsampled_times; trigger_times(i)];
        end
    end
end