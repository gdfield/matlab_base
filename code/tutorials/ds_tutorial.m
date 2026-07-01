%% DS_TUTORIAL   Identifying direction-selective (DS) retinal ganglion cells
%
% The purpose of this tutorial is to walk you through a complete, minimal
% analysis that loads a multi-electrode array (MEA) recording made during a
% drifting-grating stimulus and identifies the direction-selective (DS)
% retinal ganglion cells (RGCs) in the population.
%
% BACKGROUND
%   A direction-selective RGC fires strongly when a stimulus moves in its
%   "preferred" direction and weakly (or not at all) when the same stimulus
%   moves in the opposite, "null" direction.  To measure this we present
%   drifting gratings moving in many directions (typically 8 or 12 evenly
%   spaced angles), count the spikes each cell fires for each direction, and
%   then quantify how strongly the response depends on direction.
%
%   Two common metrics, both used below:
%     (1) Normalized vector-sum magnitude.  Treat the spike count in each
%         direction as a vector pointing in that direction; sum the vectors
%         and normalize by the total spike count.  A perfectly untuned cell
%         gives ~0 (vectors cancel); a perfectly tuned cell gives ~1.  The
%         angle of the resultant vector is the cell's PREFERRED DIRECTION.
%     (2) Direction-selectivity index (DSI) = (pref - null) / (pref + null),
%         comparing the preferred direction to the opposite direction.
%
%   DS cells stand out from the rest of the population as a cluster with high
%   values of these metrics.  We identify them by thresholding.
%
% WHAT YOU NEED
%   - A drifting-grating data set analyzed through Vision (a data0xx folder
%     containing .neurons, .params, etc.).
%   - The stimulus text file that was used to generate the gratings
%     (usually in a 'stimuli' folder next to the data, e.g. stimuli/s01.txt).
%
% This tutorial uses the lab's modern, documented DS functions:
%   get_direction_tuning, get_selectivIDX, plot_direction_tuning
% and references the classic population pipeline in the NOTES at the bottom.
%
% Created: GDF lab tutorial, 2026


%% -------------------------------------------------------------------------
%  STEP 0.  Set up the MATLAB path and the Vision Java library
%  -------------------------------------------------------------------------
%  You only need to do this once per MATLAB session.  Adjust the paths so
%  they point to wherever this code base and Vision live on your machine.

% Path to the shared lab code base (contains code/lab, etc.)
% addpath(genpath('~/Development/matlab_base/code'));

% Path to the DS/OS functions used in this tutorial
% addpath(genpath('~/Development/matlab_base/private/gfield/ds_os_functions'));

% The Vision Java library is required to read neurons/params/ei files.
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');


%% -------------------------------------------------------------------------
%  STEP 1.  Point to the data set and the stimulus file
%  -------------------------------------------------------------------------
%  EDIT THESE TWO PATHS to point at your drifting-grating recording and the
%  stimulus text file that describes it.  See ds_example.m for the pattern.

ds_data_path = '/path/to/experiment/data001/data001';      % <-- EDIT ME
stimulus_path = '/path/to/experiment/stimuli/s01.txt';     % <-- EDIT ME

% The drifting-grating stimulus is delivered in fixed-length trials.  Set the
% duration of one trial (seconds).  This is used to locate the trigger that
% marks the start of each grating presentation.  Typical values are 8 or 10 s.
trial_trigger_interval = 10;   % seconds between the start of successive trials

% Duration of the moving grating within each trial (seconds).  Spikes are
% counted over this window starting at each trigger.
grating_duration = 8;          % seconds of drifting grating per trial


%% -------------------------------------------------------------------------
%  STEP 2.  Load the data into a datarun structure
%  -------------------------------------------------------------------------
%  load_data creates the datarun structure but does not load data yet.
%  load_neurons brings in spike times (datarun.spikes) and triggers.
%  load_stim parses the stimulus file into datarun.stimulus.

datarun = load_data(ds_data_path);
datarun = load_neurons(datarun);

% Tell load_stim where the stimulus description lives, then parse it.
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_interval', trial_trigger_interval);


%% -------------------------------------------------------------------------
%  STEP 3.  Inspect the stimulus structure
%  -------------------------------------------------------------------------
%  It is worth understanding datarun.stimulus before going further:
%    .params        - the unique values of each stimulus parameter
%                     (DIRECTION, TEMPORAL_PERIOD, SPATIAL_PERIOD)
%    .combinations  - every unique combination of parameters that was shown
%    .trial_list    - which combination was shown on each trial, in order
%    .triggers      - the time (s) at which each trial began
%    .repetitions   - how many times each combination was repeated

directions = datarun.stimulus.params.DIRECTION;   % e.g. [0 45 90 ... 315]
num_dirs   = length(directions);
num_reps   = datarun.stimulus.repetitions;

fprintf('\nStimulus summary:\n');
fprintf('  %d directions: %s degrees\n', num_dirs, num2str(directions));
fprintf('  temporal periods: %s\n', num2str(datarun.stimulus.params.TEMPORAL_PERIOD));
fprintf('  spatial periods:  %s\n', num2str(datarun.stimulus.params.SPATIAL_PERIOD));
fprintf('  %d repeats per condition\n\n', num_reps);

% If the grating was run at more than one speed (temporal period), we analyze
% one speed at a time.  Here we use the first temporal/spatial period; change
% these to explore other speeds.
this_TP = datarun.stimulus.params.TEMPORAL_PERIOD(1);
this_SP = datarun.stimulus.params.SPATIAL_PERIOD(1);


%% -------------------------------------------------------------------------
%  STEP 4.  Measure each cell's direction tuning and DS metrics
%  -------------------------------------------------------------------------
%  For every cell we:
%    (a) extract the trial-by-trial spikes for each direction
%        (get_direction_tuning),
%    (b) compute the normalized vector-sum magnitude and preferred direction,
%    (c) compute the pref/null direction-selectivity index (get_selectivIDX).

cell_ids   = datarun.cell_ids;
num_cells  = length(cell_ids);

ds_magnitude   = zeros(num_cells, 1);   % normalized vector-sum magnitude [0-1]
pref_direction = zeros(num_cells, 1);   % preferred direction (degrees)
DSI            = zeros(num_cells, 1);    % (pref - null)/(pref + null)
OSI            = zeros(num_cells, 1);    % orientation-selectivity index
mean_rate      = zeros(num_cells, 1);   % mean firing rate (Hz), for QC

% Directions as a unit-complex-number per angle, used for the vector sum.
dir_unit_vectors = exp(1i * deg2rad(directions(:)'));   % 1 x num_dirs

for c = 1:num_cells

    % (a) Direction tuning for this cell at the chosen speed.
    %     spike_nums is a num_dirs x 1 vector of total spikes per direction.
    [~, spike_nums] = get_direction_tuning(datarun.spikes{c}, datarun.stimulus, ...
                                           'TP', this_TP, 'SP', this_SP, ...
                                           'stim_duration', grating_duration);
    spike_nums = spike_nums(:)';    % make it a row to match dir_unit_vectors

    total_spikes = sum(spike_nums);
    mean_rate(c) = total_spikes / (num_dirs * num_reps * grating_duration);

    % Skip essentially silent cells to avoid divide-by-zero.
    if total_spikes < eps
        continue
    end

    % (b) Normalized vector sum: magnitude tells us HOW direction selective,
    %     angle tells us the PREFERRED DIRECTION.
    vector_sum        = sum(spike_nums .* dir_unit_vectors);
    ds_magnitude(c)   = abs(vector_sum) / total_spikes;
    pref_direction(c) = mod(rad2deg(angle(vector_sum)), 360);

    % (c) Pref/null DSI (and OSI) using the lab's helper.  Works for 8 or 12
    %     evenly spaced directions.
    if num_dirs == 8 || num_dirs == 12
        [OSI(c), DSI(c)] = get_selectivIDX(spike_nums, num_dirs);
    end
end


%% -------------------------------------------------------------------------
%  STEP 5.  Visualize the population to reveal the DS cluster
%  -------------------------------------------------------------------------
%  DS cells appear as a cluster with high vector-sum magnitude and high DSI.
%  Looking at the distribution helps you choose a sensible threshold rather
%  than guessing.

figure('Name', 'DS population', 'Color', 'w');

% (1) Distribution of the vector-sum magnitude across the population.
subplot(1, 3, 1);
histogram(ds_magnitude, 0:0.02:1);
xlabel('normalized vector-sum magnitude');
ylabel('number of cells');
title('DS magnitude distribution');

% (2) Vector magnitude vs. DSI.  DS cells sit in the upper-right.
subplot(1, 3, 2);
plot(ds_magnitude, DSI, 'k.', 'MarkerSize', 10);
xlabel('vector-sum magnitude');
ylabel('DSI (pref-null)/(pref+null)');
title('magnitude vs. DSI');
axis square; box off;

% (3) Preferred directions of the most-selective cells (are they clustered
%     into the canonical DS subtypes?).  Shown as a polar scatter.
subplot(1, 3, 3);
strong = ds_magnitude > 0.3;   % preview using a provisional cutoff
polarplot(deg2rad(pref_direction(strong)), ds_magnitude(strong), 'r.', 'MarkerSize', 12);
title('preferred directions');


%% -------------------------------------------------------------------------
%  STEP 6.  Threshold to identify the DS cells
%  -------------------------------------------------------------------------
%  Choose a cutoff on the vector-sum magnitude.  Inspect the histogram in
%  Step 5 and set the threshold in the valley between the untuned bulk of the
%  population and the DS cluster.  0.3-0.5 is a common starting range.

ds_magnitude_threshold = 0.35;   % <-- TUNE ME based on the Step-5 histogram

is_ds       = ds_magnitude >= ds_magnitude_threshold;
ds_cell_ids = cell_ids(is_ds);

fprintf('Identified %d DS cells (of %d total) at magnitude threshold %.2f:\n', ...
        numel(ds_cell_ids), num_cells, ds_magnitude_threshold);
fprintf('  %s\n\n', num2str(ds_cell_ids));

% ds_cell_ids now holds the Vision cell IDs of the direction-selective RGCs.
% You can save them for downstream analysis:
% save('ds_cell_ids.mat', 'ds_cell_ids', 'ds_magnitude', 'pref_direction', 'DSI');


%% -------------------------------------------------------------------------
%  STEP 7.  Confirm by plotting example DS cells
%  -------------------------------------------------------------------------
%  Always sanity-check the identified cells.  plot_direction_tuning shows the
%  polar tuning curve (center) surrounded by spike rasters arranged by
%  direction.  A real DS cell has a clean lobe pointing toward its preferred
%  direction and near-silence at the null direction.

ds_indices = find(is_ds);
num_to_show = min(4, numel(ds_indices));

for k = 1:num_to_show
    c = ds_indices(k);

    % Re-extract the full trial-by-trial tuning (spike times, not just counts)
    % so the rasters can be drawn.
    [tuning_struct, spike_nums] = get_direction_tuning(datarun.spikes{c}, datarun.stimulus, ...
                                       'TP', this_TP, 'SP', this_SP, ...
                                       'stim_duration', grating_duration);

    fig_title = sprintf('cell %d  |  mag %.2f  |  pref %d deg', ...
                        cell_ids(c), ds_magnitude(c), round(pref_direction(c)));
    plot_direction_tuning(tuning_struct, spike_nums(:)', datarun.stimulus, ...
                          'fig_num', 10 + k, 'grating_duration', grating_duration, ...
                          'fig_title', fig_title);
end


%% -------------------------------------------------------------------------
%  NOTES AND NEXT STEPS
%  -------------------------------------------------------------------------
%
%  * MULTIPLE SPEEDS.  DS cells are often best separated when their DSI is
%    compared across two temporal periods (speeds).  The classic lab pipeline
%    does exactly this:
%        cell_ids = datarun.cell_ids;
%        [NumSpikesCell, StimComb] = get_spikescellstim(datarun, cell_ids, 0);
%        ds_struct = dscellanalysis(NumSpikesCell, StimComb);
%    ds_struct.dsindex{1} and {2} are the DSIs at the two temporal periods;
%    plotting log(dsindex{1}) vs. log(dsindex{2}) reveals the DS cluster.
%
%  * INTERACTIVE SELECTION.  DsCellFinder2.m wraps the pipeline above and lets
%    you draw a polygon around the DS cluster by hand (ginput) instead of
%    thresholding.  Useful when the cluster is not cleanly separated.
%
%  * SUBTYPES.  ON-OFF DS cells fall into ~4 preferred-direction groups and
%    ON DS cells into ~3.  Once you have ds_cell_ids, the preferred-direction
%    histogram (Step 5, panel 3) can be used to split them into subtypes.
%
%  * VALIDATION.  Cross-check DS cells against their receptive fields / EIs
%    from a companion white-noise run (see the BwPath in ds_example.m and
%    map_ei for matching cell IDs across recordings).
