function [STC, STC_diff] = compute_stc_binned(movie, spike_counts, num_frames_back, center_stimulus, batch_size)
% COMPUTE_STC_BINNED Compute spike-triggered covariance matrix from binned spike counts
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
%   STC - spike-triggered covariance matrix (n_pixels x n_pixels)
%         where n_pixels = height * width * colors * num_frames_back
%   STC_diff - STC minus stimulus covariance (excess covariance)

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
    
    % Calculate dimensionality
    if is_rgb
        n_pixels = height * width * num_colors * num_frames_back;
    else
        n_pixels = height * width * num_frames_back;
    end
    
    fprintf('Computing covariance matrix: %d x %d (%.2f MB)\n', ...
            n_pixels, n_pixels, n_pixels^2*8/(1024^2));
    
    % Compute stimulus covariance matrix
    fprintf('Computing stimulus covariance...\n');
    if is_rgb
        % Reshape movie to (n_spatial_pixels x frames)
        movie_2d = reshape(movie, height*width*num_colors, total_frames);
    else
        movie_2d = reshape(movie, height*width, total_frames);
    end
    
    % Build stimulus snippets for covariance (vectorized over valid frames)
    valid_frames = (num_frames_back+1):total_frames;
    num_valid = length(valid_frames);
    
    % Sample subset of frames for stimulus covariance to save memory
    max_stim_samples = min(num_valid, 10000);  % Use at most 10k samples
    sample_idx = randperm(num_valid, max_stim_samples);
    sampled_frames = valid_frames(sample_idx);
    
    stim_snippets = zeros(n_pixels, max_stim_samples);
    for i = 1:max_stim_samples
        frame_idx = sampled_frames(i);
        frame_indices = (frame_idx - num_frames_back + 1):frame_idx;
        if is_rgb
            snippet = movie(:, :, :, frame_indices);
            stim_snippets(:, i) = snippet(:);
        else
            snippet = movie(:, :, frame_indices);
            stim_snippets(:, i) = snippet(:);
        end
    end
    
    % Compute stimulus covariance
    stim_cov = cov(stim_snippets');
    clear stim_snippets;  % Free memory
    
    % Only process frames with spikes that have enough history
    spike_counts(1:num_frames_back) = 0;
    total_spikes = sum(spike_counts);
    
    if total_spikes == 0
        error('No valid spikes with sufficient stimulus history');
    end
    
    % Initialize accumulators for spike-triggered covariance
    sum_snippets = zeros(n_pixels, 1);
    sum_outer_products = zeros(n_pixels, n_pixels);
    
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
            
            % Vectorize snippet
            snippet_vec = snippet(:);
            
            % Accumulate mean and outer product
            sum_snippets = sum_snippets + count * snippet_vec;
            sum_outer_products = sum_outer_products + count * (snippet_vec * snippet_vec');
        end
    end
    
    % Compute spike-triggered mean
    spike_mean = sum_snippets / total_spikes;
    
    % Compute spike-triggered covariance
    % Cov = E[XX'] - E[X]E[X']
    STC = (sum_outer_products / total_spikes) - (spike_mean * spike_mean');
    
    % Compute difference from stimulus covariance
    STC_diff = STC - stim_cov;
    
    fprintf('Done!\n');
end