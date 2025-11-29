function data_list = os_adaptation_datasets(varargin)

p = inputParser;

% specify list of optional parameters
p.addParameter('experiment', 'all', @ischar);

% resolve user input and default values
p.parse(varargin{:});

path_prefix = '/Volumes/backup011/';
temp_index = 1;


%2021-09-23-0 NDF 0.0
data_list(temp_index).high_grating_datapath = [use_prefix, 'YASS/data003/data003'];
data_list(temp_index).high_stimulus_file = [use_prefix, 'stimuli/s03.txt'];
data_list(temp_index).high_trigger_interval = 10;
data_list(temp_index).low_grating_datapath = [use_prefix, 'YASS/data001/data001'];
data_list(temp_index).low_stimulus_file = [use_prefix, 'stimuli/s01.txt'];
data_list(temp_index).low_trigger_interval = 10;
data_list(temp_index).high_wn_datapath = [use_prefix, 'YASS/data002/data002'];
data_list(temp_index).low_wn_datapath = [use_prefix, 'YASS/data000/data000'];
temp_index = temp_index +1;

%2019-07-15-0
data_list(temp_index).high_grating_datapath = [path_prefix, '2019-07-15-0/data005_KR-map/data005_KR-map']; %NDF0
data_list(temp_index).high_stimulus_file = [path_prefix, '2019-07-15-0/stimuli/dg2.txt'];
data_list(temp_index).low_grating_datapath = [path_prefix, '2019-07-15-0/data003/data003']; %NDF4
data_list(temp_index).low_stimulus_path = [path_prefix, '2019-07-15-0/stimuli/dg.txt'];
data_list(temp_index).high_wn_datapath = [path_prefix, '2019-07-15-0/data004/data004']; %NDF0
data_list(temp_index).low_wn_datapath = [path_prefix, 'XXXXXXXXX'];
temp_index = temp_index + 1;


