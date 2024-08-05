function data_list = os_datasets(varargin)

p = inputParser;

% specify list of optional parameters
p.addParameter('experiment', 'all', @ischar);

% resolve user input and default values
p.parse(varargin{:});

path_prefix = '/Volumes/backup011/';
temp_index = 1;


% 2022-12-21-0 
data_list(temp_index).grating_datapath = [path_prefix, '2022-12-21-0/Yass/data001/data001'];
data_list(temp_index).stimulus_path = [path_prefix, '2022-12-21-0/stimuli/stim1.txt'];
data_list(temp_index).trigger_interval = 12;
data_list(temp_index).wn_datapath = [path_prefix, '2022-12-21-0/Yass/data000/data000'];
temp_index = temp_index +1;

% % 2023-08-09-0 
% data_list(temp_index).grating_datapath = [path_prefix, '2023-08-09-0/data001/data001'];
% data_list(temp_index).stimulus_path = [path_prefix, '2023-08-09-0 /stimuli/stim1.txt'];
% data_list(temp_index).trigger_interval = 12;
% data_list(temp_index).wn_datapath = [path_prefix, '2023-08-09-0/data000/data000'];
% temp_index = temp_index +1;


% 2012-05-31-1 
data_list(temp_index).grating_datapath = [path_prefix, '2012-05-31-1/data002/data002'];
data_list(temp_index).stimulus_path = [path_prefix, '2012-05-31-1/stimuli/s02'];
data_list(temp_index).trigger_interval = 12;
data_list(temp_index).wn_datapath = [path_prefix, '2012-05-31-1/Yass/data001/data001'];
temp_index = temp_index +1;

% 2012-10-10-1 
data_list(temp_index).grating_datapath = [path_prefix, '2012-10-10-1/data001-3600-7200/data002-map/data002-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2012-10-10-1/stimuli/s02'];
data_list(temp_index).trigger_interval = 12;
data_list(temp_index).wn_datapath = [path_prefix, '2012-10-10-1/Yass/data001-3600-7200/data001-3600-7200'];
temp_index = temp_index +1;

% 2012-10-15-0 
data_list(temp_index).grating_datapath = [path_prefix, '2012-10-15-0/data002/data002'];
data_list(temp_index).stimulus_path = [path_prefix, '2012-10-15-0/stimuli/s02'];
data_list(temp_index).trigger_interval = 12;
data_list(temp_index).wn_datapath = [path_prefix, '2012-10-15-0/Yass/data000/data000'];
temp_index = temp_index +1;

% 2012-10-31-0 
data_list(temp_index).grating_datapath = [path_prefix, '2012-10-31-0/Yass/data002/data002'];
data_list(temp_index).stimulus_path = [path_prefix, '2012-10-31-0/stimuli/s02'];
data_list(temp_index).trigger_interval = 12;
data_list(temp_index).wn_datapath = [path_prefix, '2012-10-31-0/Yass/data000/data000'];
temp_index = temp_index +1;

% %013-10-28-0
data_list(temp_index).grating_datapath = [path_prefix, '2013-10-28-0/data003/data003'];
data_list(temp_index).stimulus_path = [path_prefix, '2013-10-28-0/stimuli/S03'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2013-10-28-0/data001/data001'];
temp_index = temp_index + 1;

% 2013-10-30-0
data_list(temp_index).grating_datapath = [path_prefix, '2013-10-30-0/data003/data003'];
data_list(temp_index).stimulus_path = [path_prefix, '2013-10-30-0/stimuli/s03'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2013-10-30-0/data001/data001'];
temp_index = temp_index + 1;

% 2017-01-16-0
data_list(temp_index).grating_datapath = [path_prefix, '2017-01-16-0/data007/data007'];
data_list(temp_index).stimulus_path = [path_prefix, '2017-01-16-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2017-01-16-0/data006_KR/data006_KR'];
temp_index = temp_index + 1;

% 2017-06-07-0
data_list(temp_index).grating_datapath = [path_prefix, '2017-06-07-0/data007-map/data007-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2017-06-07-0/stimuli/data007.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2017-06-07-0/data005/data005'];
temp_index = temp_index + 1;

%2017-07-20-0 DG and WN
data_list(temp_index).grating_datapath = [path_prefix, '2017-07-20-0/data007-map/data007-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2017-07-20-0/stimuli/data007.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2017-07-20-0/data004_KR/data004_KR'];
temp_index = temp_index + 1;

% %2017-09-08-0 DG and WN
data_list(temp_index).grating_datapath = [path_prefix, '2017-09-08-0/data006-map/data006-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2017-09-08-0/stimuli/data006.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2017-09-08-0/data005/data005'];
temp_index = temp_index + 1;

%2018-02-26-0 DG and WN
data_list(temp_index).grating_datapath = [path_prefix, '2018-02-26-0/data004-map/data004-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2018-02-26-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-02-26-0/data003/data003'];
temp_index = temp_index + 1;

% %2018-04-09-0 DG and WN
% data_list(temp_index).grating_datapath = [path_prefix, '2018-04-09-0/data004-map/data004-map'];
% data_list(temp_index).stimulus_path = [path_prefix, '2018-04-09-0/stimuli/dg.txt'];
% data_list(temp_index).trigger_interval = 10;
% data_list(temp_index).wn_datapath = [path_prefix, '2018-04-09-0/data005/data005'];
% temp_index = temp_index + 1;
% 

%2018-04-18-0 DG and WN
data_list(temp_index).grating_datapath = [path_prefix, '2018-04-18-0/data003-map_KR/data003'];
data_list(temp_index).stimulus_path = [path_prefix, '2018-04-18-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-04-18-0/data002_KR/data002_KR'];
temp_index = temp_index + 1;

%2018-05-23-0 DG and WN
data_list(temp_index).grating_datapath = [path_prefix, '2018-05-23-0/data003-map/data003-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2018-05-23-0/stimuli/data003.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-05-23-0/data002/data002'];
temp_index = temp_index + 1;

% %2018-09-10-0
data_list(temp_index).grating_datapath = [path_prefix, '2018-09-10-0/data006-map/data006-map'];
data_list(temp_index).stimulus_path = [path_prefix, '2018-09-10-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-09-10-0/data005/data005']; 
temp_index = temp_index + 1;

%2018-11-26-0
data_list(temp_index).grating_datapath = [path_prefix, '2018-11-26-0/data004-map/data004-map']; %NDF0
data_list(temp_index).stimulus_path = [path_prefix, '2018-11-26-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-11-26-0/data003/data003']; %NDF0
temp_index = temp_index + 1;

% %2018-11-29-0
data_list(temp_index).grating_datapath = [path_prefix, '2018-11-29-0/data013-map/data013-map']; %NDF0
data_list(temp_index).stimulus_path = [path_prefix, '2018-11-29-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-11-29-0/data010_KR/data010_KR']; %NDF0, this WN run has not been analyzed
temp_index = temp_index + 1;

%% potentially version good. Worthly of YASS or Kilosort work.
% %2018-11-30-0
data_list(temp_index).grating_datapath = [path_prefix, '2018-11-30-0/data004_KR-map/data004_KR-map']; %NDF0
data_list(temp_index).stimulus_path = [path_prefix, '2018-11-30-0/stimuli/dg.txt'];
data_list(temp_index).trigger_interval = 10;
data_list(temp_index).wn_datapath = [path_prefix, '2018-11-30-0/data003/data003']; %NDF0
temp_index = temp_index + 1;








% 
% 
% 
% 
% % 2024-01-17-0
% grating_datapath = [path_prefix, 'Rat/2024-01-17-0/data001/'];
% stimulus_path = [path_prefix, 'Rat/2024-01-17-0/stimuli/stim1.txt'];
% wn_datapath = [path_prefix, 'Rat/2024-01-17-0/data000'];
% 
% %2019-07-15-0
% grating_datapath = [path_prefix,'2019-07-15-0/data005_KR-map/data005_KR-map']; %NDF0
% stimulus_path = [path_prefix,'2019-07-15-0/stimuli/dg2.txt'];
% %grating_datapath = [path_prefix,'2019-07-15-0/data003/data003']; %NDF4
% %stimulus_path = [path_prefix,'2019-07-15-0/stimuli/dg.txt'];
% wn_datapath = [path_prefix,'2019-07-15-0/data004/data004']; %NDF0
% 
% 



% %% SC rats
% 
% % 2021-07-13-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-0/data000/data000';
% 
% % 2021-07-13-1
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-13-1/data000/data000';
% 
% % 2021-07-15-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-07-15-0/data000/data000';
% 
% % 2021-06-09-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/data002/data002';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/stimuli/s02.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-09-0/data000/data000';
% 
% % 2021-06-02-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/data002/data002';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/stimuli/s02.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-06-02-0/data001/data001';
% 
% %% OS experiments
% 
% %2021-09-09-0 NDF 3.0 VISION
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data000/data000';
% 
% %2021-09-09-0 NDF 3.0 YASS
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data000/data000';
% 
% %2021-09-09-0 NDF 0.0 VISION
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s03.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data002/data002';
% 
% %2021-09-09-0 NDF 0.0 YASS
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s03.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data002/data002';
% 
% %2021-09-09-0 NDF 0.0, GABA block
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data004/data004';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s04.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data005/data005';
% 
% %2021-09-09-0 NDF 0.0, GABA block, YASS
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data004/data004';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s04.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/YASS/data005/data005';
% 
% %2021-09-09-0 NDF, GABA block washout
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data006/data006';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/stimuli/s06.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-09-0/data007/data007';
% 
% %2021-09-23-0 NDF 3.0 VISION
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data000/data000';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s00.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data002/data002';
% 
% %2021-09-23-0 NDF 0.0 VISION
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s03.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/data004/data004';
% 
% %2021-09-23-0 NDF 3.0 YASS
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/YASS/data000/data000';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s00.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/YASS/data002/data002';
% 
% %2021-09-23-0 NDF 0.0 YASS
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/YASS/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/stimuli/s03.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-09-23-0/YASS/data004/data004';
% 
% %2021-10-07-0 Yass NDF0.0 control
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-10-07-0/Yass/data001/data001';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-10-07-0/stimuli/s01.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-10-07-0/Yass/data000/data000'
% 
% %% Others from notebooks -- light adaptation
% 
% %2019-07-15-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/data005_KR-map/data005_KR-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/stimuli/dg2.txt';
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/data003/data003'; %NDF4
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/data004/data004'; %NDF0
% 
% %2019-04-01-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-01-0/data002_MR-map/data002_MR-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-01-0/stimuli/dg2.txt'; %change stimulus path and run 
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-01-0/data001/data001'; %NDF0
% 
% %2019-02-15-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-02-15-0/data011-map/data011-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-02-15-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-02-15-0/data009-map/data009-map'; %NDF0
% 
% %2019-01-21-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/data004-map/data004-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/data003/data003'; %NDF0
% 
% %2018-11-30-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/data004_KR-map/data004_KR-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/data003/data003'; %NDF0
% 
% %2018-11-29-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-29-0/data013-map/data013-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-29-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-29-0/data012/data012'; %NDF0, this WN run has not been analyzed
% 
% %2018-11-26-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/data004-map/data004-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/data003/data003'; %NDF0
% 
% %2018-09-10-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-09-10-0/data006-map/data006-map'; %NDF0
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-09-10-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-09-10-0/data005/data005'; %NDF0
% 
% %% hi density array data sets
% 
% %2013-10-28-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-28-0/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-28-0/stimuli/S03';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-28-0/data001/data001';
% 
% %2013-10-30-0
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-30-0/data003/data003';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-30-0/stimuli/s03';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-30-0/data001/data001';
% 
% %% Datasets from Kiersten's files
% 
% %2018-02-26-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-02-26-0/data004-map/data004-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-02-26-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-02-26-0/data003/data003';
% 
% %2018-04-09-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-09-0/data004-map/data004-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-09-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-09-0/data005/data005';
% 
% %2018-04-18-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-18-0/data003-map_KR/data003-map_KR';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-18-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-18-0/data002_KR/data002_KR';
% 
% %2018-05-23-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-05-23-0/data003-map/data003-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-05-23-0/stimuli/data003.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-05-23-0/data002/data002';
% 
% %2018-11-26-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/data004-map/data004-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-26-0/data003/data003';
% 
% %2018-11-30-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/data004_KR-map/data004_KR-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-11-30-0/data003_KR/data003_KR';
% 
% %2019-01-21-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/data004_KR-map/data004_KR-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-01-21-0/data003_KR/data003_KR';
% 
% %2019-04-08-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-08-0/data002_KR-map/data002_KR-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-08-0/stimuli/dg.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-04-08-0/data001/data001';
% 
% %2019-07-15-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/data005_KR-map/data005_KR-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/stimuli/dg2.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-15-0/data004_KR/data004_KR';
% 
% %% contrast reversing gratings
% 
% %2017-07-20-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-07-20-0/data007-map/data007-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-07-20-0/stimuli/data007.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-07-20-0/data004_KR/data004_KR';
% 
% %2017-09-08-0 DG and WN
% grating_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-09-08-0/data006-map/data006-map';
% stimulus_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-09-08-0/stimuli/data006.txt';
% wn_datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2017-09-08-0/data005/data005';


