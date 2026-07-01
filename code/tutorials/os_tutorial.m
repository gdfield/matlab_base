%% OS_TUTORIAL   Identifying orientation-selective (OS) retinal ganglion cells
%
% This tutorial is the companion to ds_tutorial.m.  It loads a drifting-
% grating multi-electrode array (MEA) recording and identifies the
% ORIENTATION-selective (OS) retinal ganglion cells (RGCs) in the population.
%
% BACKGROUND
%   An orientation-selective RGC responds strongly to a grating (or bar) of a
%   particular ORIENTATION and weakly to the orthogonal orientation.  Unlike a
%   direction-selective (DS) cell, a pure OS cell does not care which of the
%   two opposite directions the grating drifts: a vertical grating drifting
%   left and the same vertical grating drifting right have the same
%   orientation, so they drive an OS cell equally.
%
%   This 180-degree periodicity is the crucial difference from DS analysis:
%
%     DIRECTION is periodic over 360 deg  -> DS uses a vector sum with the
%                                            actual angle theta:      exp(i*theta)
%     ORIENTATION is periodic over 180 deg -> OS uses a "double-angle"
%                                            vector sum:              exp(i*2*theta)
%
%   Doubling the angle maps directions 180 deg apart onto the SAME angle, so
%   the two lobes of an orientation-tuned response add constructively instead
%   of cancelling.  We therefore:
%     (1) fold the direction tuning into orientation tuning (average the two
%         opposite directions), and
%     (2) compute a normalized double-angle vector-sum magnitude.  The angle
%         of the resultant, halved, is the cell's PREFERRED ORIENTATION.
%
%   Metric summary (both computed below):
%     - Normalized double-angle vector-sum magnitude: 0 (untuned) to 1 (sharp).
%     - Orientation-selectivity index (OSI) = (pref - orth)/(pref + orth),
%       comparing the preferred orientation to the orthogonal orientation
%       (get_selectivIDX).
%
%   IMPORTANT CONFOUND.  Because gratings drift, a DS cell also produces an
%   apparently oriented response (it fires for one direction, which lies on
%   one orientation axis).  So DS cells will show up with a high OS magnitude
%   too.  To isolate cells that are genuinely orientation- but NOT direction-
%   selective, we exclude cells with a high DS magnitude (Step 6).
%
% This tutorial uses the lab's modern DS/OS functions:
%   get_direction_tuning, get_selectivIDX, plot_direction_tuning
%
% Created: GDF lab tutorial, 2026


%% -------------------------------------------------------------------------
%  STEP 0.  Set up the MATLAB path and the Vision Java library
%  -------------------------------------------------------------------------
%  You only need to do this once per MATLAB session.

% addpath(genpath('~/Development/matlab_base/code'));
% addpath(genpath('~/Development/matlab_base/private/gfield/ds_os_functions'));
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');


%% -------------------------------------------------------------------------
%  STEP 1.  Point to the data set and the stimulus file
%  -------------------------------------------------------------------------
%  EDIT THESE TWO PATHS to point at your drifting-grating recording and the
%  stimulus text file that describes it.

ds_data_path  = '/path/to/experiment/data001/data001';     % <-- EDIT ME
stimulus_path = '/path/to/experiment/stimuli/s01.txt';     % <-- EDIT ME

trial_trigger_interval = 10;   % seconds between the start of successive trials
grating_duration       = 8;    % seconds of drifting grating per trial


%% -------------------------------------------------------------------------
%  STEP 2.  Load the data into a datarun structure
%  -------------------------------------------------------------------------

datarun = load_data(ds_data_path);
datarun = load_neurons(datarun);

datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_interval', trial_trigger_interval);


%% -------------------------------------------------------------------------
%  STEP 3.  Inspect the stimulus structure
%  -------------------------------------------------------------------------

directions = datarun.stimulus.params.DIRECTION;   % e.g. [0 45 90 ... 315]
num_dirs   = length(directions);
num_reps   = datarun.stimulus.repetitions;

fprintf('\nStimulus summary:\n');
fprintf('  %d directions: %s degrees\n', num_dirs, num2str(directions));
fprintf('  temporal periods: %s\n', num2str(datarun.stimulus.params.TEMPORAL_PERIOD));
fprintf('  spatial periods:  %s\n', num2str(datarun.stimulus.params.SPATIAL_PERIOD));
fprintf('  %d repeats per condition\n\n', num_reps);

% Analyze one speed at a time; use the first temporal/spatial period.
this_TP = datarun.stimulus.params.TEMPORAL_PERIOD(1);
this_SP = datarun.stimulus.params.SPATIAL_PERIOD(1);

% OS analysis assumes directions come in opposite pairs (each direction has a
% partner 180 deg away), which holds for the usual 8- or 12-direction runs.
if mod(num_dirs, 2) ~= 0
    warning('OS analysis expects an even number of directions in opposite pairs.');
end


%% -------------------------------------------------------------------------
%  STEP 4.  Measure each cell's orientation tuning and OS metrics
%  -------------------------------------------------------------------------
%  For every cell we:
%    (a) extract spikes per direction (get_direction_tuning),
%    (b) compute the normalized DOUBLE-ANGLE vector-sum magnitude and the
%        preferred orientation,
%    (c) compute the OSI (get_selectivIDX),
%    (d) also compute the DS magnitude so we can exclude DS cells later.

cell_ids  = datarun.cell_ids;
num_cells = length(cell_ids);

os_magnitude     = zeros(num_cells, 1);   % normalized double-angle magnitude [0-1]
pref_orientation = zeros(num_cells, 1);   % preferred orientation (0-180 deg)
OSI              = zeros(num_cells, 1);    % (pref - orth)/(pref + orth)
ds_magnitude     = zeros(num_cells, 1);   % single-angle magnitude, to reject DS cells
mean_rate        = zeros(num_cells, 1);   % mean firing rate (Hz), for QC

% Unit vectors for the two conventions.
dir_unit_vectors    = exp(1i * deg2rad(directions(:)'));       % single angle  -> DS
orient_unit_vectors = exp(1i * 2 * deg2rad(directions(:)'));   % double angle  -> OS

for c = 1:num_cells

    % (a) Direction tuning for this cell at the chosen speed.
    [~, spike_nums] = get_direction_tuning(datarun.spikes{c}, datarun.stimulus, ...
                                           'TP', this_TP, 'SP', this_SP, ...
                                           'stim_duration', grating_duration);
    spike_nums = spike_nums(:)';    % row vector, aligned with directions

    total_spikes = sum(spike_nums);
    mean_rate(c) = total_spikes / (num_dirs * num_reps * grating_duration);

    if total_spikes < eps
        continue
    end

    % (b) Orientation: double-angle vector sum.  The magnitude measures how
    %     orientation selective the cell is; half the resultant angle is the
    %     preferred orientation (wrapped into 0-180 deg).
    orient_vector_sum   = sum(spike_nums .* orient_unit_vectors);
    os_magnitude(c)     = abs(orient_vector_sum) / total_spikes;
    pref_orientation(c) = mod(rad2deg(angle(orient_vector_sum)) / 2, 180);

    % (c) OSI via the lab helper (also folds opposite directions internally).
    if num_dirs == 8 || num_dirs == 12
        OSI(c) = get_selectivIDX(spike_nums, num_dirs);
    end

    % (d) DS magnitude (single-angle) so we can separate OS from DS below.
    ds_magnitude(c) = abs(sum(spike_nums .* dir_unit_vectors)) / total_spikes;
end


%% -------------------------------------------------------------------------
%  STEP 5.  Visualize the population to reveal the OS cluster
%  -------------------------------------------------------------------------

figure('Name', 'OS population', 'Color', 'w');

% (1) Distribution of the orientation (double-angle) vector-sum magnitude.
subplot(1, 3, 1);
histogram(os_magnitude, 0:0.02:1);
xlabel('normalized orientation magnitude');
ylabel('number of cells');
title('OS magnitude distribution');

% (2) OS magnitude vs. DS magnitude.  This is the key panel for separating
%     the two: pure OS cells have HIGH OS magnitude but LOW DS magnitude
%     (upper-left); DS cells drift toward the diagonal / lower-right.
subplot(1, 3, 2);
plot(ds_magnitude, os_magnitude, 'k.', 'MarkerSize', 10);
xlabel('DS magnitude (single-angle)');
ylabel('OS magnitude (double-angle)');
title('OS vs. DS');
axis square; box off;

% (3) Preferred orientations of the strongly-tuned cells.  OS RGC types tend
%     to cluster around cardinal (0/90) or oblique orientations.  Plotted on a
%     0-180 deg axis (doubled so both lobes are visible on the polar plot).
subplot(1, 3, 3);
strong = os_magnitude > 0.25;   % provisional cutoff, just for this preview
polarplot(deg2rad(2 * pref_orientation(strong)), os_magnitude(strong), ...
          'b.', 'MarkerSize', 12);
title('preferred orientations (angle doubled)');


%% -------------------------------------------------------------------------
%  STEP 6.  Threshold to identify the OS cells (excluding DS cells)
%  -------------------------------------------------------------------------
%  Choose a cutoff on the orientation magnitude from the Step-5 histogram, and
%  ALSO require the DS magnitude to be low so we keep cells that are
%  orientation- but not direction-selective.

os_magnitude_threshold = 0.30;   % <-- TUNE ME from the Step-5 histogram
ds_magnitude_maximum   = 0.20;   % <-- reject cells that are mostly DS

is_os       = (os_magnitude >= os_magnitude_threshold) & ...
              (ds_magnitude <= ds_magnitude_maximum);
os_cell_ids = cell_ids(is_os);

fprintf('Identified %d OS cells (of %d total)\n', numel(os_cell_ids), num_cells);
fprintf('  OS magnitude >= %.2f and DS magnitude <= %.2f:\n', ...
        os_magnitude_threshold, ds_magnitude_maximum);
fprintf('  %s\n\n', num2str(os_cell_ids));

% To include ALL orientation-tuned cells (DS cells included), simply drop the
% DS constraint:
%   is_os = os_magnitude >= os_magnitude_threshold;

% Save for downstream analysis:
% save('os_cell_ids.mat', 'os_cell_ids', 'os_magnitude', 'pref_orientation', 'OSI');


%% -------------------------------------------------------------------------
%  STEP 7.  Confirm by plotting example OS cells
%  -------------------------------------------------------------------------
%  A genuine OS cell shows TWO lobes 180 deg apart in the polar tuning plot
%  (strong response to both drift directions along the preferred orientation)
%  and near-silence at the orthogonal orientation.  Contrast this with a DS
%  cell, which shows a single lobe.

os_indices  = find(is_os);
num_to_show = min(4, numel(os_indices));

for k = 1:num_to_show
    c = os_indices(k);

    [tuning_struct, spike_nums] = get_direction_tuning(datarun.spikes{c}, datarun.stimulus, ...
                                       'TP', this_TP, 'SP', this_SP, ...
                                       'stim_duration', grating_duration);

    fig_title = sprintf('cell %d  |  OS mag %.2f  |  pref ori %d deg', ...
                        cell_ids(c), os_magnitude(c), round(pref_orientation(c)));
    plot_direction_tuning(tuning_struct, spike_nums(:)', datarun.stimulus, ...
                          'fig_num', 20 + k, 'grating_duration', grating_duration, ...
                          'fig_title', fig_title, 'color', 'b');
end


%% -------------------------------------------------------------------------
%  NOTES AND NEXT STEPS
%  -------------------------------------------------------------------------
%
%  * OS vs. DS.  The single most useful diagnostic is the OS-vs-DS scatter
%    (Step 5, panel 2).  Pure OS cells populate the high-OS / low-DS corner.
%    If you instead want every orientation-tuned cell regardless of direction
%    tuning, drop the ds_magnitude constraint in Step 6.
%
%  * WHY DOUBLE THE ANGLE.  Orientation has a period of 180 deg, so a cell
%    tuned to "vertical" fires for gratings drifting both up and down.  Using
%    exp(i*2*theta) folds those two directions onto the same phase so their
%    responses reinforce; the resultant's angle, halved, is the preferred
%    orientation.  (This is the circular-statistics standard for axial data.)
%
%  * SHARPNESS / TUNING WIDTH.  Beyond a binary OS/not-OS call, you can fit the
%    folded orientation tuning curve with a von Mises function to estimate
%    tuning width; see ComputeMises.m in the dscellanalysis/ folder.
%
%  * VALIDATION.  Cross-check OS cells against receptive fields / EIs from a
%    companion white-noise run (see ds_example.m BwPath and map_ei) and, if
%    available, confirm the preferred-orientation clustering expected for the
%    OS RGC types in your species.
