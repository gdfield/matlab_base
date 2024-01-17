function params=wholeArrayAutoAnalysisLargeScale(eiFilePath, pathToPreparation,folders,pathSave,varargin)
% Uses Gonzalo's spike sorting code (GM2) to automatically analyze all the
% neurons in a preparation.
% Inputs:       eiFilePath (e.g. '/Volumes/Analysis/2015-10-06-3/data000/data000.ei')
%               pathToPreparation (e.g. '/Volumes/Analysis/2015-10-06-3/data001-data002/')
%               folders (e.g.  {'data001','data002','data003'})
%   optional    arrayID: for now is an optional parameter but we could make
%                   it a default param. Options
%               threshold
%               sortingThreshold
%               startPattern
%               endPattern
% GM 04/2016 (BASED ON wholeArrayAutoAnalysis LG 12/2015)


% Set the default values for optional parameters
arrayID = 501;
threshold = 10; % Analyze any cell whose EI is greater than this value on the stimulating electrode
sortingThreshold = 45; % Analyzed cells must have an EI of at least this magnitude on the recording electode
patternList=[1:512];
override = 0;
params0=[];
sortSmallCell=0;
electrodes=[];
Art0=[];
useKernels=1;
cellIDs0=[];
useBrownian = 1;
templates0=[];
updateVar=0;
% Read the optional input arguments
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'patternlist'
            patternList = varargin{j*2};
        case 'threshold'
            threshold = varargin{j*2};
        case 'sortingthreshold'
            sortingThreshold = varargin{j*2};
        case 'arrayid'
            arrayID = varargin{j*2};
        case 'params'
            params0= varargin{j*2};
        case 'sortsmallcell'
            sortSmallCell=varargin{j*2};
        case 'cellids'
            cellIDs0=varargin{j*2};
        case 'usekernels'
            useKernels= varargin{j*2};
        case 'templates'
            templates0 = varargin{j*2};
        case 'usebrownian'
            useBrownian =varargin{j*2};
        case 'initialart'
            Art0= varargin{j*2};
            case 'updatevar'
            updateVar= varargin{j*2};
        case 'electrodes'
            electrodes=varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

arrayObj = Array(arrayID);
numElecs = arrayObj.getNumElecs;

if numElecs ~= 512 && override ==0
    endPattern = 519;
end

% Initialize the analysis
if(isempty(params0))
    params = InitializeArrayRobust(pathToPreparation,arrayObj,'foldernames',folders);
    
    % Adjust parameters?
    params.global.sortData=1;
    params.global.nTrial=80;
    params.global.subSampleRate=1;
    params.bundle.findBundle=1;
    params.bundle.useBundleAlg=1;
    params.global.useStimElec=0;
    params.global.useStimElectrodeBeforeBundle=1;
    params.global.useStimElectrodeAfterBundle=0;
    params.global.thresEI=40;
    params.global.saveArt=1;
    params.path.pathToEi=eiFilePath;
else
    params=params0;
    clear params0
end
params.global.sortSmallCell=sortSmallCell;

% Find pattern-cell pairs to analyze.
disp('Finding pattern-cell pairs to analyze...');
if(isempty(cellIDs0))
[eiM,nIdList] = convertEiFileToMatrix(eiFilePath);
cellsToAnalyze = cell(numElecs,1);
end
listElectrodes=patternList;
if(~isempty(electrodes))
    listElectrodes=electrodes;
end

for e = listElectrodes
    electrodeNo = e;
    if(isempty(cellIDs0))
    cellIDs = [];
    for n = 1:1:size(nIdList,1)
        tempID  = nIdList(n);          % gets the id number of the ith neuron
        ei = squeeze(eiM(:,:,n))';
        eiMin = abs(min(ei));
        if eiMin(electrodeNo) > threshold && (max(eiMin(:)) > sortingThreshold)
            cellIDs = cat(1,cellIDs,tempID);
        end
    end
    if(isempty(cellIDs) && sortSmallCell)
        maximums=max(abs(squeeze(eiM(electrodeNo,:,:))));

    cellIDs=nIdList(find(maximums==max(maximums)));
    end
    else
        cellIDs =cellIDs0;
    end
    cellsToAnalyze{e} = cellIDs;
    fprintf('*');
    if mod(e,80) == 0
        fprintf('electrode %0.0f \n',e);
    end
    
end

if(~isempty(electrodes))
    cells=[];
    for e=electrodes
        cells=cat(1,cells,cellsToAnalyze{e});
    end
    clear cellsToAnalyze
    cellsToAnalyze{patternList(1)}=unique(cells);
end
fprintf('\n');

% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%
% %% %% % Calculate the total number of pattern-cell pairs to analyze % %%
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%
totalNo = 0;
for e = patternList
    totalNo = length(cellsToAnalyze{e}) + totalNo;
end
fprintf('Total number of pattern-cell pairs is %0.0f\n',totalNo);

% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%
% %% %% % Update to reflect the latest pattern-cell pair
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%


% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%
% %% %% % Begin auto analysis % %% %% %
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %%
disp('Beginning automatic analysis. Come back in a while');


for e = patternList
    
    cellIDs = cellsToAnalyze{e} ;
    patternNo = e;
    disp( '**************************************************');
    disp(['****** analyzing pattern no. ' num2str(e) ' *****************']);
    disp( '**************************************************');
    if ~isempty(cellIDs)
        if(isempty(cellIDs0))
        templates=makeTemplatesFromEiShift(eiFilePath, cellIDs,arrayObj,[1:numElecs]);
        else
            templates=templates0;
        end
        try
            elecStimSortingLargeScale(pathToPreparation,folders,pathSave,e,cellIDs,templates,params,arrayObj,'usekernels',useKernels,'initialart',Art0,'useBrownian',useBrownian,'updateVar',updateVar)
        catch 
            disp(['Unable to analyze pattern ' num2str(patternNo)])
            %rethrow(ME);
            
        end
    end
    
    try
        disp(['finished analyzing pattern '...
            num2str(patternNo) ' neurons ' num2str(cellIDs') '. pattern ' num2str(e-patternList(1)+1) ' of ' num2str(endPattern-startPattern+1)']);
    catch
        try
            disp(['finished analyzing pattern '...
                num2str(patternNo) ' neurons ' num2str(cellIDs) '. pattern ' num2str(e-patternList(1)+1) ' of ' num2str(endPattern-startPattern+1)']);
        catch
            disp(['finished analyzing pattern ' num2str(patternNo) ]);
        end
    end
end

%catch ME
%   disp(['Final electrode analyzed was ' num2str(e)])
%  rethrow(ME)
%end