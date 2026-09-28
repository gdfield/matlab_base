function datarun = os_trigger_check(datarun, stimulus_path, trigger_interval) 
%
% usage: datarun = os_trigger_check(datarun, stimulus_path) 
%
% This function handles some user issues in setting up the triggers for
% drifting gratings during experiments in the Field Lab under the Photons
% stimulation code.
%
% The function takes a datarun structure and hands back a datarun structure
% with a stimulus field that contains the trigger times and other
% information associated with a drifting grating run. Some datasets have an
% unusual trigger logic (b/c of user error) that need to be handled
% case-by-case. This function, takes care of that.
%
% Created: GDF 2026-01-08
%

if contains(datarun.names.rrs_prefix, '2021-09-09-0/data001/')
    % handle these fucked up triggers
    set_trig_one = find(diff(datarun.triggers) > 1.8 & diff(datarun.triggers) < 2.3);
    set_trig_two = find(diff(datarun.triggers) > 4);
    set_trig_three = find(diff(datarun.triggers) > 1.01 & diff(datarun.triggers) < 1.1);
    all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
    all_trig_indices_sorted = sort(all_trig_indices, 'ascend');

    % load stimulus information for gratings
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);

elseif contains(datarun.names.rrs_prefix, '2021-09-23-0/data003')
    % handle these fucked up triggers
    all_trig_indices_sorted = [1, find(diff(datarun.triggers) > 2.01)'];
    
    % load stimulus information for gratings
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);

elseif contains(datarun.names.rrs_prefix, '2022-12-21-0/Yass/data001')
    for i=1:length(datarun.triggers)-1
        trig_dif(i) = round(datarun.triggers(i+1)-datarun.triggers(i));
    end
    trial_trig = find(trig_dif==2);
    trial_trig = [1,trial_trig+1];
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim_amr(datarun,'user_defined_trigger_set', trial_trig);
elseif contains(datarun.names.rrs_prefix, '2024-08-06-0/data002') || ...
       contains(datarun.names.rrs_prefix, '2024-08-13-0/data002')
    % Recording started ~2.6 s into trial 1 (7 of 9 cycle triggers present),
    % and the trial 1->2 boundary gap is ~1.4 s, so the 2-s-gap rule merged
    % trials 1 and 2 and shifted every later trial by one stimulus condition.
    % Fix: detect boundaries as any non-cycle interval, require every trial
    % after trial 1 to match the stimulus file, then drop the first
    % repetition (trial 1 is truncated; get_direction_tuning needs equal reps).
    % 2026-09-27 GDF + Claude
    [trial_trig, rep] = find_amr_trial_onsets(datarun.triggers, stimulus_path);
    if rep.n_onsets ~= rep.n_file_trials || ~isequal(rep.count_mismatch_trials, 1)
        error('os_trigger_check: onset validation failed for %s (onsets %d, file trials %d)', ...
            datarun.names.rrs_prefix, rep.n_onsets, rep.n_file_trials);
    end
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim_amr(datarun, 'user_defined_trigger_set', trial_trig);
    nc = length(datarun.stimulus.combinations);
    keep = nc+1:length(datarun.stimulus.trials);
    datarun.stimulus.trials      = datarun.stimulus.trials(keep);
    datarun.stimulus.trial_list  = datarun.stimulus.trial_list(keep);
    datarun.stimulus.triggers    = datarun.stimulus.triggers(keep);
    datarun.stimulus.repetitions = datarun.stimulus.repetitions - 1;
    fprintf('os_trigger_check: validated %d trials vs stim file; dropped truncated first repetition -> %d reps\n', ...
        rep.n_onsets, datarun.stimulus.repetitions);
    
elseif contains(datarun.names.rrs_prefix, '2024-11-14-0/data001')
    % One trigger was dropped inside trial 161 (a TP-60 trial), leaving a
    % ~2 s within-trial gap that the 2-s-gap rule read as a trial boundary.
    % That split one trial into two (385 markers for 384 trials) and shifted
    % every trial after 161 by one stimulus condition. find_amr_trial_onsets
    % removes the spurious split; the true onset of trial 161 is unaffected,
    % so all 8 repetitions are kept.
    % 2026-09-27 GDF + Claude
    [trial_trig, rep] = find_amr_trial_onsets(datarun.triggers, stimulus_path);
    mm = rep.count_mismatch_trials;
    if rep.n_onsets ~= rep.n_file_trials || ~isequal(mm, rep.merged_trials) || ...
            any(rep.expected_counts(mm) - rep.observed_counts(mm) ~= 1) || ...
            ~isempty(rep.irregular_intervals)
        error('os_trigger_check: onset validation failed for %s', datarun.names.rrs_prefix);
    end
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim_amr(datarun, 'user_defined_trigger_set', trial_trig);
    fprintf('os_trigger_check: validated %d trials vs stim file; repaired dropped trigger in trial(s) %s\n', ...
        rep.n_onsets, mat2str(rep.merged_trials));
    
elseif contains(datarun.names.rrs_prefix, '2024-01-17-0/data001')
    % Triggers were clean all along (432 regular trials); the stimulus file in
    % this dataset's folder was the wrong one (a copy of 2023-11-08-0's).
    % new_os_datasets.m now points to the matching file. Require a perfect
    % per-trial match to the stimulus file before loading.
    % 2026-09-27 GDF + Claude
    [trial_trig, rep] = find_amr_trial_onsets(datarun.triggers, stimulus_path);
    if rep.n_onsets ~= rep.n_file_trials || ~isempty(rep.count_mismatch_trials) || ...
            ~isempty(rep.merged_trials) || ~isempty(rep.irregular_intervals)
        error('os_trigger_check: onset validation failed for %s', datarun.names.rrs_prefix);
    end
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim_amr(datarun, 'user_defined_trigger_set', trial_trig);
    fprintf('os_trigger_check: validated %d trials vs stim file (exact match)\n', rep.n_onsets);
    
elseif any([...
       contains(datarun.names.rrs_prefix, '2023-07-19-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-07-27-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-10-04-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-10-20-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-10-27-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-10-31-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2023-11-08-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2024-01-30-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2024-02-19-0/data001'), ...
       contains(datarun.names.rrs_prefix, '2024-05-08-0/data001')])
    % 2023-2024 datasets on the newer rig protocol (540-frame stimuli), using
    % the 2-second-gap trial-marker rule from 2022-12-21-0 above.
    % STATUS 2026-09-27: NOT fully validated. Total trial counts match for most,
    % but a count match does not guarantee correct trial->condition pairing.
    % (2024-01-17-0 moved to its own branch: wrong stimulus file, see new_os_datasets.m.)
    % Candidates for migration to find_amr_trial_onsets (see branch above).
    for i=1:length(datarun.triggers)-1
        trig_dif(i) = round(datarun.triggers(i+1)-datarun.triggers(i));
    end
    trial_trig = find(trig_dif==2);
    trial_trig = [1,trial_trig+1];
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim_amr(datarun,'user_defined_trigger_set', trial_trig);

else
    % load stimulus information for gratings
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim(datarun, 'user_defined_trigger_interval', trigger_interval);
end
