function [STV, STV_diff] = compute_stv_binned(movie, spike_counts, num_frames_back, center_stimulus, batch_size)
% COMPUTE_STV_BINNED Compute spike-triggered variance from binned spike counts
%
% Inputs:
%   movie - 4D array (height x width x colors x frames) for RGB
%           or 3D array (height x width x frames) for grayscale
%   spike_counts - vector of spike counts per frame (length = frames)
%   num_frames_back - number of frames to include before each spike
%   center_stimulus - (optional) if true, subtract mean from movie (default: true)
%   batch_size - (optional) number of spike frames to process per batch (default: inf)
%
% Outputs:
%   STV - spike-triggered variance (height x width x colors x num_frames_back)
%         or (height x width x num_frames_back) for grayscale
%   STV_diff - STV minus stimulus variance (excess variance)

    if nargin < 4
        center_stimulus = false;
    end
    if nargin < 5
        batch_size = inf;  % Process all spikes at once by default
    end

    % Handle both RGB (4D) and grayscale (3D) inputs
    movie_dims = size(movie);
    if length(movie_dims) == 4
        % RGB movie: height x width x colors x frames
        [height, width, num_colors, total_frames] = size(movie);
        is_rgb = true;
    elseif length(movie_dims) == 3
        % Grayscale movie: height x width x frames
        [height, width, total_frames] = size(movie);
        num_colors = 1;
        is_rgb = false;
    else
        error('Movie must be 3D (grayscale) or 4D (RGB)');
    end
    
    % Mean-center the stimulus if requested
    if center_stimulus
        movie = movie - mean(movie(:));
    end
    
    % Verify spike_counts matches movie length
    if length(spike_counts) ~= total_frames
        error('spike_counts must have same length as movie frames');
    end
    
    % Compute stimulus variance before modifying spike_counts
    if is_rgb
        stim_var = var(movie, 0, 4);  % variance along 4th dimension (frames)
    else
        stim_var = var(movie, 0, 3);  % variance along 3rd dimension (frames)
    end
    
    % Only process frames with spikes that have enough history
    spike_counts(1:num_frames_back) = 0;
    total_spikes = sum(spike_counts);
    
    if total_spikes == 0
        error('No valid spikes with sufficient stimulus history');
    end
    
    % Initialize accumulators
    if is_rgb
        sum_snippets = zeros(height, width, num_colors, num_frames_back);
        sum_sq_snippets = zeros(height, width, num_colors, num_frames_back);
    else
        sum_snippets = zeros(height, width, num_frames_back);
        sum_sq_snippets = zeros(height, width, num_frames_back);
    end
    
    % Find frames with spikes
    spike_frames = find(spike_counts > 0);
    num_spike_frames = length(spike_frames);
    
    % Process in batches to avoid memory issues
    num_batches = ceil(num_spike_frames / batch_size);
    
    if num_batches > 1
        fprintf('Processing %d spikes in %d batches...\n', total_spikes, num_batches);
    end
    
    for batch = 1:num_batches
        % Define batch range
        start_idx = (batch - 1) * batch_size + 1;
        end_idx = min(batch * batch_size, num_spike_frames);
        batch_spike_frames = spike_frames(start_idx:end_idx);
        
        if num_batches > 1
            fprintf('  Batch %d/%d: processing frames %d to %d\n', ...
                    batch, num_batches, start_idx, end_idx);
        end
        
        % Accumulate for this batch
        for i = 1:length(batch_spike_frames)
            spike_frame = batch_spike_frames(i);
            count = spike_counts(spike_frame);
            frame_indices = (spike_frame - num_frames_back + 1):spike_frame;
            
            if is_rgb
                snippet = movie(:, :, :, frame_indices);
            else
                snippet = movie(:, :, frame_indices);
            end
            
            sum_snippets = sum_snippets + count * snippet;
            sum_sq_snippets = sum_sq_snippets + count * (snippet.^2);
        end
    end
    
    % Compute STA (needed for variance calculation)
    STA = sum_snippets / total_spikes;
    
    % Compute STV and STV_diff
    STV = (sum_sq_snippets / total_spikes) - STA.^2;
    STV_diff = STV - stim_var;
end