%OS experiments
%2021-09-09-0 NDF 3.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data000/data000';


%2021-09-09-0 NDF 3.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data001/data001';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data000/data000';

%2021-09-09-0 NDF 0.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data003/data003';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s03.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data002/data002';




%2021-09-09-0 NDF 0.0
grating_datapath = '/Users/gdfield/Analysis/2021-09-23-0/data003/data003';
stimulus_path = '/Users/gdfield/Analysis/2021-09-23-0//stimuli/s03.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data002/data002';




% load stuff
datarun = load_data(grating_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');


set_trig_one = find(diff(datarun.triggers) > 1.8 & diff(datarun.triggers) < 2.3);
set_trig_two = find(diff(datarun.triggers) > 4);
set_trig_three = find(diff(datarun.triggers) > 1.01 & diff(datarun.triggers) < 1.1);
all_trig_indices = [1, set_trig_one', set_trig_two', set_trig_three'];
all_trig_indices_sorted = sort(all_trig_indices, 'ascend');

% load stimulus information for gratings
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);



%%
%2021-09-23-0 NDF 3.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data000/data000';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s00.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data001/data001';


% load stuff
datarun = load_data(grating_datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');

all_trig_indices_sorted = [1, find(diff(datarun.triggers) > 2.01)'];

% load stimulus information for gratings
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_set', all_trig_indices_sorted);




%%


%2021-09-23-0 NDF 3.0
grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data003/data003';
stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s03.txt';
wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data004/data004';


