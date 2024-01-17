function elecStimSortingLargeScale(pathToPreparation,dataFolders,pathSave,patternNos,neuronIds,templates,params,arrayObj,varargin)
% inputs:     pathToAnalysisData: '/Volumes/Analysis/...'choose the electrical stim data folder (001,002, etc) that contains data organized by pattern
%             pathToEi: a string, '/full/path/to/ei/file.ei'
%             patternNos
%             neuronIds

useKernels=1;
Art0=[];
useBrownian=1;
updateVar=0;
% codebase_path = matlab_code_path;
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
       
         case 'usekernels'
            useKernels= varargin{j*2};
        case 'initialart'
            Art0= varargin{j*2};
        case 'usebrownian'
            useBrownian = varargin{j*2};
        case 'updatevar'
            updateVar =varargin{j*2};
    
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end


for p = 1:length(patternNos)
    patternNo=patternNos(p);   
    [elecRespAuto]=DoSpikeSortingLargeScale(pathToPreparation,patternNo,neuronIds,params,templates,arrayObj,'FoldersNames',dataFolders,'usekernels',useKernels,'initialart',Art0,'useBrownian',useBrownian);
    fname = fullfile(pathSave,['elecRespAuto_p' ...
    num2str(elecRespAuto.stimInfo.patternNo) '.mat']);   
    %num2str(elecRespAuto.neuronInfo.neuronIds') '_p' ...
       
    %fname = fullfile(pathSave,['elecRespAuto_n' ...
        %num2str(elecRespAuto.neuronInfo.neuronIds') '_p' ...
        
    save(fname,'elecRespAuto');
    disp(['done analyzing ' fname]);
end




