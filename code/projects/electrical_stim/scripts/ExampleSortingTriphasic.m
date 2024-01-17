%% Input parameters (preparation, etc)
pathToPreparation='/Volumes/Analysis/2012-09-24-3/'; %set preparation folder
folders={'data003','data004','data005','data006'}; %set folder(s) to look for stimulation patterns (px folders)

eiFilePath='/Volumes/Analysis/2012-09-24-3/data000/data000.ei'; %directory to look for templates
pathSave=[pathToPreparation 'Autosort-single/']; %The folders where responses will be saved
mkdir(pathSave) %create diretory to save data
arrayObj=Array(501); %set 512 array
patternList=[1:512];

%% Initialize model
[params]=InitializeArrayRobust(pathToPreparation,arrayObj,'foldernames',folders,'useBrownian',0,'useLocalization',[1 1]);

% Algorithm parameters. Don't touch this (they will change so that they
% produce optimal results)
params.bundle.findBundle=1;
params.global.sortData=1;
params.global.nTrial=80;
params.global.subSampleRate=1;
params.bundle.useBundleAlg=0;
params.global.useStimElectrodeBeforeBundle=1;
params.global.useStimElectrodeAfterBundle=0;
params.global.filterStimElectrode=1;
params.global.extraStimElectrode=1;
params.global.filterNoStimElectrodes=1;
params.global.extraNoStimElectrodes=0;
params.global.saveArt=0;


%% Do the actual spike sorting
wholeArrayAutoAnalysisLargeScale(eiFilePath,pathToPreparation,folders,pathSave,'patternList',patternList,'arrayid',501,'useBrownian',0,'updateVar',0,'params',params);
