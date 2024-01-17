%% Input parameters (preparation, etc)
pathToPreparation='/Volumes/Analysis/2015-04-09-2/'; %set preparation folder
folders={'data003'}; %set folder(s) to look for stimulation patterns (px folders)

eiFilePath='/Volumes/Analysis/2015-04-09-2/data001/data001.ei'; %directory to look for templates
pathSave=[pathToPreparation 'Autosort-single/']; %The folders where responses will be saved
mkdir(pathSave) %create diretory to save data
arrayObj=Array(501); %set 512 array

patternList=[ 1 2 15 16 39 40 47 48 55 56 57 58 64 67 68 71 72 75 76 77 78 93 94 109 110 135 136 137 138 139 140 143 144 159 160 161 162 163 178 179 180
    

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
params.global.useStimElectrodeAfterBundle=1;
params.global.filterStimElectrode=0;
params.global.extraStimElectrode=0;
params.global.filterNoStimElectrodes=0;
params.global.extraNoStimElectrodes=0;
params.global.saveArt=0;


%% Do the actual spike sorting
wholeArrayAutoAnalysisLargeScale(eiFilePath,pathToPreparation,folders,pathSave,'patternList',patternList,'arrayid',501,'useBrownian',0,'updateVar',0,'params',params,'useKernels',0);
