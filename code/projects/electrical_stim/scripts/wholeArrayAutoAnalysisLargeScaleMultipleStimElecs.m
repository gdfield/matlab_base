function wholeArrayAutoAnalysisLargeScaleMultipleStimElecs(eiFilePath, pathToPreparation,folder,pathLoad,pathSave,startPattern,endPattern,params,varargin)
% Given a folder with elecRespAuto files corresponding to single-electrode
% stimulation, does multi-electrode stimulation spikeSorting
% 
% GM 07/2016 (BASED ON wholeArrayAutoAnalysis LG 12/2015)

arrayID = 501; 
threshold = 10; % Analyze any cell whose EI is greater than this value on the stimulating electrode 
sortingThreshold = 45; % Analyzed cells must have an EI of at least this magnitude on the recording electode
override = 0; 
params0=[];
sortSmallCell=0;
electrodes=[];
Art0=[];
useKernels=1;
cellIDs0=[];
templates0=[];
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
       case 'arrayid'
            arrayID = varargin{j*2};
         case 'cellids'
            cellIDs0=varargin{j*2};
        case 'usekernels'
            useKernels= varargin{j*2};
        case 'templates'
            templates0 = varargin{j*2};
            
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end


dirs=dir(pathLoad);

for p=startPattern:endPattern
try
    patternNo=p;
    pathToAnalysisData=[pathToPreparation folder{1}];
[listAmpsNew listStimElecs TracesAll4 Art4 channelsConnected]=loadAmps(pathToAnalysisData,patternNo);
[listAmpsNew2 listAmpsNew] = collapseAmplitudesEqual(listAmpsNew, listAmpsNew);


stimElecs=listStimElecs(1,:);

clear  Art0, 
clear fArt0
clear amps
clear xStims
clear KersStims
clear polarity
        
for e=1:length(stimElecs)
    
    str=['p' num2str(stimElecs(e)) '.mat'];
    l=length(str);
    for j=1:length(dirs)
        try
            
        if(isequal(dirs(j).name(end-l+1:end),str))
           
        load([pathLoad dirs(j).name])
        Art0{e}=elecRespAuto.stimInfo.Art;
        fArt0{e}=elecRespAuto.stimInfo.firstArt;
        amps{e}=elecRespAuto.stimInfo.listAmps';
        xStims{e}=elecRespAuto.stimInfo.xStim;
        KersStims{e}=elecRespAuto.stimInfo.KersStim;
        polarity(e)=elecRespAuto.stimInfo.polarity;
        
        
        end
    end
    end
end
    
arrayObj=Array(arrayID);    
clear Apred
for e=1:length(stimElecs);
stimElec=stimElecs(e);
Kers{1}=params.patternInfo.Kers{1};
Kers{3}=params.patternInfo.Kers{3};

Diag=arrayObj.difPos(setdiff(arrayObj.getElectrodes(),stimElec),setdiff(arrayObj.getElectrodes(),stimElec))/arrayObj.maxR;
positions=arrayObj.getPositions;
[~,rho] = cart2pol(positions(:,1)-positions(stimElec,1),positions(:,2)-positions(stimElec,2));
rho=rho(setdiff(arrayObj.getElectrodes(),stimElec));
[Ker, KerD]=evalKernels(Diag,rho/arrayObj.maxR,params.arrayInfo.x(4:6),1);
Kers{2}=Ker;

[Apred1]=ExtrapolateArtifactForMultiple(Art0{e}(:,:,1:params.global.Tmax),params.patternInfo.var0,params.arrayInfo.x,Kers,setdiff(arrayObj.getElectrodes(),stimElecs(e)),abs(amps{e}),polarity(e),abs(listAmpsNew(:,e)),sign(listAmpsNew(1,e)),max(params.arrayInfo.listAmpsAll));
[Apred2]=ExtrapolateArtifactForMultipleStimElec(Art0{e}(:,:,1:params.global.Tmax),params.patternInfo.var0,xStims{e},KersStims{e},abs(amps{e}),polarity(e),abs(listAmpsNew(:,e)),sign(listAmpsNew(1,e)),max(params.arrayInfo.listAmpsAll),stimElec);

Apred(:,stimElec,:)=Apred2;
Apred(:,setdiff(arrayObj.getElectrodes(),stimElec),:)=Apred1;

Art00{e}=Apred;

end


ArtInitial=0;
for e=1:length(stimElecs)
    ArtInitial=ArtInitial+Art00{e};
end
for e=1:length(stimElecs)
    ArtInitial(1,stimElecs(e),:)=ArtInitial(1,stimElecs(e),:)+reshape(squeeze(fArt0{e}(1,1,stimElecs(e),:)),1,1,params.global.Tmax);
end

params.global.saveArt=0;
wholeArrayAutoAnalysisLargeScale(eiFilePath, pathToPreparation,folder,pathSave,'startpattern',patternNo,'endpattern',patternNo,'arrayid',arrayID,'params',params,'usekernels',0,'initialart',ArtInitial,'electrodes',stimElecs,'templates',templates0,'cellids',cellIDs0);

end


end



end


