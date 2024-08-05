%% Check OSI and DSI for ReaChR cells 

%% load relevant file paths

addpath(genpath('/Users/marijarudzite/Desktop/matlab/code/lab/')); %path to matlab/code/lab folder 
addpath(genpath('/Users/marijarudzite/Desktop/matlab/code/projects/')); %path to matlab/code/projects folder 
addpath(genpath('/Users/marijarudzite/Documents/MATLAB/DS&OS stuff/ds_os_functions/')); %path to ds code stuff
addpath(genpath('//Users/marijarudzite/Documents/MATLAB/Retinotectal_analysis_for_AMR/'));%path to folder with ReaChR analysis specific functions 

% display format 
format shortg; 
%% data and save paths 
reach_data_path = '/Users/marijarudzite/Desktop/SC_stuff/RAT Data/';

rat = 'R23-07';
experiment = '2023-10-04-0'; %folder name of the experiment 
YASS = false; %set to true for YASS false for Vison

% Specify data filenames
dfiles = 'data001/data001'; 
d_folder = 'data001';

%data path
datapath = strcat(reach_data_path,rat,'/',experiment,'/'); % path to analysis folder on server 
stimulus_path = strcat(datapath,'stimuli/stim1.txt');
if YASS
    dpath = strcat(datapath,'Yass/');
else %vision
    dpath = datapath;
end

%save path 
matlab_path = strcat(reach_data_path,rat,'/Matlab_analysis/');

%% load cells to plot
load(strcat(matlab_path,'mapped_cells.mat'));

%% load data 

ei_params = struct('array_id',551,'array_type',512); % {id, type} = {1551; 519} [dense] or, {551; 512} [sparse]
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei',1, 'load_ei_params', ei_params);
dataruns = load_data(strcat(dpath,dfiles), opt);
dataruns.piece.array_type = ei_params.array_type;
dataruns.piece.array_id = ei_params.array_id;
[dataruns,array_info] = load_array_info(dataruns);


%% load stimulus
%mannually find trial star triggers
for i=1:length(dataruns.triggers)-1
trig_dif(i) = round(dataruns.triggers(i+1)-dataruns.triggers(i));
end
trial_trig = find(trig_dif==2);
trial_trig = [1,trial_trig+1];

%load stimulus 
dataruns.names.stimulus_path = stimulus_path;
dataruns = load_stim_amr(dataruns,'user_defined_trigger_set', trial_trig)

% extract spatial and temporal periods from datarun
spatial_periods = dataruns.stimulus.params.SPATIAL_PERIOD;
temp_periods = dataruns.stimulus.params.TEMPORAL_PERIOD;
num_sps = length(spatial_periods);
num_tps = length(temp_periods);
num_dirs = length(dataruns.stimulus.params.DIRECTION);
