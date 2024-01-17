%% Input parameters (preparation, etc)
pathToPreparation='/Volumes/Analysis/2016-04-21-6/'; %set preparation folder
folders={'data003'}; %set folder(s) to look for stimulation patterns (px folders)

eiFilePath='/Volumes/Analysis/2016-04-21-6/data001/data001.ei'; %directory to look for templates
pathSave=[pathToPreparation 'Autosort-single/']; %The folders where responses will be saved
mkdir(pathSave) %create diretory to save data
arrayObj=Array(501); %set 512 array
patternList=[508 482 437 384 363 176 347 327 328 89 282 190 181 190 112 69];

%% Initialize model
[params]=InitializeArrayRobust(pathToPreparation,arrayObj,'foldernames',folders,'useBrownian',0,'useLocalization',[1 1]);

% Algorithm parameters. Don't touch this (they will change so that they
% produce optimal results)
params.bundle.findBundle=1;
params.global.sortData=1;
params.global.nTrial=80;
params.global.subSampleRate=1;
params.bundle.useBundleAlg=0;
params.global.useStimElectrodeBeforeBundle=0;
params.global.useStimElectrodeAfterBundle=0;
params.global.filterStimElectrode=0;
params.global.extraStimElectrode=0;
params.global.filterNoStimElectrodes=1;
params.global.extraNoStimElectrodes=1;
params.global.saveArt=0;
params.global.extraStimElectrode=0;

%% Do the actual spike sorting
wholeArrayAutoAnalysisLargeScale(eiFilePath,pathToPreparation,folders,pathSave,'patternList',patternList,'arrayid',501,'useBrownian',0,'updateVar',0,'params',params);
