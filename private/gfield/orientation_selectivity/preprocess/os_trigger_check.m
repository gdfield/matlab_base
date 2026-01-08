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
else
    % load stimulus information for gratings
    datarun.names.stimulus_path = stimulus_path;
    datarun = load_stim(datarun, 'user_defined_trigger_interval', trigger_interval);
end