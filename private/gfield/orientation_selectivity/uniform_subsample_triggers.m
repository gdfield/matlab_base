function uniform_times = uniform_subsample_triggers(trigger_times, samples_per_interval)
% UNIFORM_SUBSAMPLE_TRIGGERS Create uniformly spaced samples between triggers
%
% Inputs:
%   trigger_times - vector of trigger times (frame indices)
%   samples_per_interval - number of samples per interval between triggers
%
% Outputs:
%   uniform_times - vector of uniformly spaced sample times
%                   Length = (num_triggers - 1) * samples_per_interval

    % Sort trigger times
    trigger_times = sort(trigger_times(:));
    
    if length(trigger_times) < 2
        error('Need at least 2 trigger times to create intervals');
    end
    
    num_intervals = length(trigger_times) - 1;
    uniform_times = zeros(num_intervals * samples_per_interval, 1);
    
    idx = 1;
    
    % For each pair of consecutive triggers
    for i = 1:num_intervals
        start_time = trigger_times(i);
        end_time = trigger_times(i + 1);
        
        % Create samples_per_interval equally spaced points in this interval
        % linspace with samples_per_interval points between (but not including) endpoints
        interval_samples = linspace(start_time, end_time, samples_per_interval + 2);
        interval_samples = interval_samples(2:end-1);  % Exclude endpoints
        
        % Store in output array
        uniform_times(idx:idx+samples_per_interval-1) = interval_samples;
        idx = idx + samples_per_interval;
    end
end