function data_list = new_os_datasets()
% new_os_datasets   Dataset registry for 20 of the 22 originally-identified
% unanalyzed datasets. 2023-08-24-0 and 2023-10-17-0 were excluded: neither
% has a multi-direction drifting-grating stimulus on disk (their assigned
% files are either single-direction SP/contrast sweeps or a different
% stimulus class entirely).
% Returns the same struct format as os_datasets() -- fields:
%   grating_datapath, stimulus_path, trigger_interval, wn_datapath
%
% Paths marked REVIEW should be verified against the disk layout before
% running curation. Once curation is complete and entries are confirmed,
% paste them into os_datasets.m manually.
%
% 2026-09-26 generated for os_cell_finder_script_claude_auto.m

path_prefix = '/Volumes/gdf/rat-data/';
temp_index = 1;

% ==========================================================
% Clear pairings (scanner auto-matched dataNNN to stimulus)
% ==========================================================

% 2018-07-18-0
data_list(temp_index).grating_datapath = [path_prefix, '2018-07-18-0/data008-map/data008-map'];
data_list(temp_index).stimulus_path    = [path_prefix, '2018-07-18-0/stimuli/s08.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2018-07-18-0/data007-map/data007-map'];
temp_index = temp_index + 1;

% 2020-03-19-0  % dg009 confirmed as the 12-direction OS grating (dg000-dg008 are single-direction SP/contrast sweeps); data009-map has full sorted output
data_list(temp_index).grating_datapath = [path_prefix, '2020-03-19-0/data009-map/data009-map'];
data_list(temp_index).stimulus_path    = [path_prefix, '2020-03-19-0/stimuli/dg009.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2020-03-19-0/data006/data006'];
temp_index = temp_index + 1;

% 2021-07-13-0
data_list(temp_index).grating_datapath = [path_prefix, '2021-07-13-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2021-07-13-0/stimuli/s01.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2021-07-13-0/data000/data000'];
temp_index = temp_index + 1;

% 2021-07-13-1
data_list(temp_index).grating_datapath = [path_prefix, '2021-07-13-1/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2021-07-13-1/stimuli/s01.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2021-07-13-1/data000/data000'];
temp_index = temp_index + 1;

% 2021-07-15-0
data_list(temp_index).grating_datapath = [path_prefix, '2021-07-15-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2021-07-15-0/stimuli/s01.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2021-07-15-0/data000/data000'];
temp_index = temp_index + 1;

% 2021-12-28-0  % REVIEW: also has s03/s04/s06; confirm s01/data001 is the drifting-grating run
data_list(temp_index).grating_datapath = [path_prefix, '2021-12-28-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2021-12-28-0/stimuli/s01.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2021-12-28-0/data000/data000'];
temp_index = temp_index + 1;

% 2024-02-19-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-02-19-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-02-19-0/stimuli/s01.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-02-19-0/data000/data000'];
temp_index = temp_index + 1;

% ==========================================================
% Inferred pairings (stim1/stim2 naming -- confirm dataNNN <-> stimN)
% ==========================================================

% 2023-07-19-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-07-19-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-07-19-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-07-19-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-07-27-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-07-27-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-07-27-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-07-27-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-10-04-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-10-04-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-10-04-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-10-04-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-10-20-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-10-20-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-10-20-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-10-20-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-10-27-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-10-27-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-10-27-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-10-27-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-10-31-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-10-31-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-10-31-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-10-31-0/data000/data000'];
temp_index = temp_index + 1;

% 2023-11-08-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2023-11-08-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2023-11-08-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2023-11-08-0/data000/data000'];
temp_index = temp_index + 1;

% 2024-01-17-0  | NOTE: the stimuli/ folder here holds a copy of 2023-11-08-0's files (generator stim_2023-11-08-0.m; stim1.txt byte-identical to 2023-11-08-0/stim1.txt) and does NOT match this recording. The recorded trial sequence matches 2024-02-19-0/stimuli/s01.txt (identical to 2023-10-31-0/stim1.txt): 432/432 per-trial trigger counts, and within-condition reliability 0.41 vs 0.02 for the local file. Verified 2026-09-27 GDF + Claude.
data_list(temp_index).grating_datapath = [path_prefix, '2024-01-17-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-02-19-0/stimuli/s01.txt'];  % see note above
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-01-17-0/data000/data000'];
temp_index = temp_index + 1;

% 2024-01-30-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-01-30-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-01-30-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-01-30-0/data000/data000'];
temp_index = temp_index + 1;

% 2024-05-08-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-05-08-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-05-08-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-05-08-0/data000/data000'];
temp_index = temp_index + 1;

% 2024-08-06-0  % REVIEW: stimuli has stim2/stim3 only (no stim1)  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-08-06-0/data002/data002'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-08-06-0/stimuli/stim2.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-08-06-0/data001/data001'];
temp_index = temp_index + 1;

% 2024-08-13-0  % REVIEW: stim2/stim3 pattern like 2024-08-06-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-08-13-0/data002/data002'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-08-13-0/stimuli/stim2.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-08-13-0/data001/data001'];
temp_index = temp_index + 1;

% 2024-11-14-0  | NOTE: 540-frame/newer-rig protocol -- loaded via os_trigger_check's load_stim_amr branch; trigger_interval below is unused for this dataset
data_list(temp_index).grating_datapath = [path_prefix, '2024-11-14-0/data001/data001'];
data_list(temp_index).stimulus_path    = [path_prefix, '2024-11-14-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath      = [path_prefix, '2024-11-14-0/data000/data000'];
temp_index = temp_index + 1;

% ==========================================================
% DEFERRED: 2023-08-09-0 -- no stimuli/ folder on disk.
% ==========================================================
