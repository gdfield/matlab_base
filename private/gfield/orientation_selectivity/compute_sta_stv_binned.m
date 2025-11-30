function [STA, STV] = compute_sta_stv_binned(movie, spike_counts, num_frames_back, center_stimulus)
% COMPUTE_STA_STV_BINNED Compute STA/STV from binned spike counts
%
% Inputs:
%   movie - 4D array (height x width x colors x frames) for RGB
%           or 3D array (height x width x frames) for grayscale
%   spike_counts - vector of spike counts per frame (length = frames)
%   num_frames_back - number of frames to include before each spike
%   center_stimulus - (optional) if true, subtract mean from movie (default: true)
%
% Outputs:
%   STA - spike-triggered average (height x width x colors x num_frames_back)
%         or (height x width x num_frames_back) for grayscale
%   STV - spike-triggered variance (height x width x colors x num_frames_back)
%         or (height x width x num_frames_back) for grayscale

    if nargin < 4
        center_stimulus = true;
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
    
    % Accumulate weighted by spike count
    for i = 1:length(spike_frames)
        spike_frame = spike_frames(i);
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
    
    % Compute STA and STV
    STA = sum_snippets / total_spikes;
    STV = (sum_sq_snippets / total_spikes) - STA.^2;
end