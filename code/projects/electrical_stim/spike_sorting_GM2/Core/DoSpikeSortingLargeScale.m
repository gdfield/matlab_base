function [Output]=DoSpikeSortingLargeScaleTest(pathToPreparation,patternNo,neuronIds,params,templates,arrayObj,varargin)
% Gonzalo Mena, 03/2016
% Optional arguments:
%                     'FoldersNames': default, empty and then searches all
%                     folders for desired pattern. Otherwise specify a cell
%                     array of folders to search through to find the
%                     pattern
% Lauren Grosberg 6/2016 edits to handle the array object
% Set optional arguments.
FoldersNames = [];
useKernels=1;
TracesSimulated=[];
sampledTrials=[];
meanTraces=[];
onsetInput=[];
pvals=[];
Art0=[];
updateVar=0;
useBrownian=0;
useBrownianStimElec=0;
nsubAmplitudes=[];
% Read in optional inputs.
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
        case 'foldersnames'
            FoldersNames = varargin{j*2};
        case 'usekernels'
            useKernels= varargin{j*2};
        case 'initialart'
            Art0= varargin{j*2};
        case 'pvals'
            pvals=varargin{j*2};
        case 'simulation'
            TracesSimulated= varargin{j*2};
        case 'usebrownian'
            useBrownian= varargin{j*2};
        case 'usebrownianstimelec'
            useBrownianStimElec= varargin{j*2};
        case 'meantraces'
            meanTraces=varargin{j*2};;
        case 'onsetinput'
            onsetInput=varargin{j*2};
        case 'sampledtrials'
            sampledTrials=varargin{j*2};
        case 'updatevar'
            updateVar=varargin{j*2};
        case 'xpostbundle'
            xPostBundle=varargin{j*2};
        case 'nsubamplitudes'
            nsubAmplitudes=varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

%positions=params.global.positions;
x=params.arrayInfo.x;
useBundleAlg = params.bundle.useBundleAlg;
useStimElectrodeBeforeBundle = params.global.useStimElectrodeBeforeBundle;
useStimElectrodeAfterBundle = params.global.useStimElectrodeAfterBundle;

Tmax=params.global.Tmax;
nTrial=params.global.nTrial;

if isempty(FoldersNames)
    dirs=dir(pathToPreparation);
    cont=1;
    for i=1:length(dirs)
        if(length(dirs(i).name)>=4)
            aux=find(dirs(i).name(1:4)=='data',1);
            if(~isempty(aux))
                FoldersNames{cont}=dirs(i).name; %#ok<AGROW>
                cont=cont+1;
            end
        end
    end
end

for f=1:length(FoldersNames)
    pathAux=[pathToPreparation FoldersNames{f}];
    dirs=dir(pathAux);
    for i=1:length(dirs)
        aux=find(dirs(i).name(1)=='p',1);
        if(~isempty(aux))
            if(isequal(dirs(i).name(2:end),num2str(patternNo)))
                path=pathAux;
            end
        end
    end
end

pathToAnalysisData=[path filesep];



[TracesAll, Art, var01, listAmps, listCurrents, stimElecs, onset, onsetC, pval, Res, sampledTrials, Trace0, polarity firstArt channelsConnected]=loadTracesArtSort(pathToAnalysisData,patternNo,Tmax,nTrial,params.global.subSampleRate,arrayObj,'params',params,'sampledtrials',sampledTrials,'nsubamplitudes',nsubAmplitudes);

var0=params.patternInfo.var0;
if(~isempty(onsetInput))
    onsetC=onsetInput;
    pval=pvals;
    if(~isnan(onsetC))
        onset=listAmps(onsetC);
    else
        onset=NaN;
    end
end

if(~isempty(TracesSimulated))
    clear Art
    TracesAll=TracesSimulated{1};
    listAmps=TracesSimulated{2};
    listCurrents=TracesSimulated{3};
    for i=1:size(TracesAll,1)
        Art(i,:,:)=squeeze(nanmean(TracesAll(i,:,:,:),2));
    end
end

%stimElecs=1; %to make it work with p2

%shift stimulating elecrode to match the new ordering inducing after
%erasing inactive channels




% Get active electrode list from the array object
activeElecs = arrayObj.getElectrodes;
numActiveElecs = length(activeElecs);
positions=arrayObj.getPositions;
[~,rho] = cart2pol(positions(:,1)-positions(stimElecs(1),1),positions(:,2)-positions(stimElecs(1),2));



ind=setdiff(activeElecs,find(rho==0));
if(length(stimElecs)>1)
    ind=setdiff(activeElecs,stimElecs);
end


for e=1:length(stimElecs)
    stimElecsRel(e)=find(activeElecs==stimElecs(e));
end
indRel=setdiff([1:numActiveElecs],stimElecsRel);
params.patternInfo.ind=ind;
params.patternInfo.indRel=indRel;
params.patternInfo.stimElecs=stimElecs;
params.patternInfo.stimElecsRel=stimElecsRel;


if size(Art,2) ~= numActiveElecs
    disp('size of the artifact should match the number of active channels');
end



if(~isempty(meanTraces))
    for i=1:size(TracesAll,1)
        
        TracesAll(i,:,:,:)=TracesAll(i,:,:,:)-repmat(nanmean(TracesAll(i,:,:,:),2),1,size(TracesAll,2),1,1);
        
    end
end




if(~isempty(Art0))
    for i=1:size(TracesAll,1)
        try
            TracesAll(i,:,stimElecsRel,:)=TracesAll(i,:,stimElecsRel,:)+repmat(firstArt(i,stimElecsRel,:),size(TracesAll,2),1,1);
            
            TracesAll(i,:,:,:)=TracesAll(i,:,:,:)-reshape(repmat(Art0(i,:,:),size(TracesAll,2),1,1),1,size(TracesAll,2),size(TracesAll,3),size(TracesAll,4));
            
        end
    end
end



Dif1= zeros(Tmax);
for i=1:Tmax
    for j=1:Tmax
        Dif1(i,j)=abs(i-j);
    end
end

if(useBrownian==0)
    Dif3 = zeros(size(listAmps,1),size(listAmps,1));
    for j=1:length(listAmps)
        for i=1:length(listAmps)
            Dif3(i,j)=abs(listAmps(i)-listAmps(j));
        end
    end
else
    for j=1:length(listAmps)
        for i=1:length(listAmps)
            Dif3(i,j)=min(listAmps(i),listAmps(j));
        end
    end
end
rho=rho(ind);

Difs{2}=arrayObj.difPos(ind,ind)/arrayObj.maxR;
Difs{3}=Dif3;
Difs{3}=Difs{3}/max(max(params.arrayInfo.listAmpsAll));
Dif1=Dif1/Tmax;
Difs{1}=Dif1;


Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/arrayObj.maxR;
Diags{3}=listAmps(:,1)/max(params.arrayInfo.listAmpsAll);

DifsWhiten=Difs;
DiagsWhiten=Diags;
DifsWhiten{2}=params.arrayInfo.Dif2(activeElecs,activeElecs)/arrayObj.maxR;
DiagsWhiten{2}=ones(numActiveElecs,1);
DiagsWhiten{1}=ones(Tmax,1);

var0old=var0;

if(params.global.whiten)
    [TracesWhiten templatesWhiten var0 varsE]=WhitenTracesTemplates(TracesAll,templates,params,DifsWhiten,DiagsWhiten,[]);
else
    TracesWhiten=TracesAll;
    templateswhiten=templates;
end

if(params.global.whiten)
    
    for i=1:length(templatesWhiten)
        for e=1:length(activeElecs)
            elsWhiten(i,e)=max(abs(templatesWhiten{i}(e,:)));
            els(i,e)= max(abs(templates{i}(e,:)));
        end
    end
    
    [a b c d e]=regress(elsWhiten(:),[ones(length(elsWhiten(:)),1), els(:)]);
    
    params.global.thresEI=30*a(2);
    templates=templatesWhiten;
    TracesAll=TracesWhiten;
else
    params.global.thresEI=30;
end

if(useKernels)
    
    [KersWhitenCholeskyInv]=WhitenKernels(Difs,Diags,params.arrayInfo.xwhiten,size(Art,1));
end


%change this later!!!!
if(isnan(onsetC))
    onsetC=size(TracesAll,1)+1;
end



if(useBrownian)
    types=[1 1 8];
    factp=[0 3 6 7];
else
    types=[1 1 1];
    factp=[0 3 6 9];
    
end

if(useKernels==1&&updateVar==1)
    if(onsetC>1)
        DifAux=Difs;
        DifAux{3}=Dif3(1:onsetC-1,1:onsetC-1);
        DiagsAux=Diags;
        DiagsAux{3}=listAmps(1:onsetC-1)'/max(params.arrayInfo.listAmpsAll);
        if(useBrownian)
            f1=@(Art,y)logDetKron(Art(1:onsetC-1,ind,:),[x(1:7) y(1) log(var0old)],DifAux,8,types,DiagsAux,[3 3 1]);
            g1=@(y)f1(Art,y);
            x1 = fminunc(g1,[x(end)],params.global.options);
            x([8])=x1;
        else
            f1=@(Art,y)logDetKron(Art(1:onsetC-1,ind,:),[x(1:9) y(1) log(var0old)],DifAux,[10],types,DiagsAux,[3 3 3]);
            g1=@(y)f1(Art,y);
            x1 = fminunc(g1,[x(end)],params.global.options);
            x([10])=x1;
        end
        params.arrayInfo.x=x;
    end
end


if(useKernels)
    for k=1:3
        [Ker, KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),types(k));
        if(k==2);
            if(  params.global.whiten)
                
                Ker=diag(1./sqrt(varsE(indRel)))*Ker*diag(1./sqrt(varsE(indRel)));
            end
        end
        Kers{k}=KersWhitenCholeskyInv{k}*Ker*KersWhitenCholeskyInv{k}';
        [a, b]=eig(Kers{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
        
    end
    
    x(end)=x(end)-params.arrayInfo.xwhiten(3);
    params.arrayInfo.x=x;
    
end

if(params.global.whiten)
    [Art templatesWhiten]=WhitenTracesTemplates(reshape(Art,1,size(Art,1),size(Art,2),size(Art,3)),[],params,DifsWhiten,DiagsWhiten,varsE);
    Art=squeeze(Art);
    ArtWhiten=Art;
else
    ArtWhiten=Art;
end


%this should not be declared here but at initializearrayRobust
params.global.useBrownian=useBrownian;

% Update params structure.
disp(['Update params structue for pattern ' num2str(patternNo)]);
if(useKernels)
    params.patternInfo.Q=Q;
    params.patternInfo.Qt=Qt;
    params.patternInfo.Kers=Kers;
    params.patternInfo.dL=dL;
end
params.patternInfo.Art =ArtWhiten;
params.patternInfo.var0=var0;
params.patternInfo.rho=rho;
params.patternInfo.patternNo=patternNo;
params.patternInfo.listAmps=listAmps;
params.patternInfo.listCurrents=listCurrents;
params.patternInfo.Diags=Diags;
params.patternInfo.Difs=Difs;
params.bundle.onsBundle=onsetC;
params.bundle.onset=onset;

params.neuronInfo.templates = templates;
breakpoints=findBreakStimElecs(listCurrents);
params.patternInfo.breakpoints{1}=[0 breakpoints{1}' size(Art,1)];
        
if(params.global.filterStimElectrode+params.global.extraStimElectrode>=1)
    if(useKernels)
        
        Diff = zeros(size(listAmps,1),size(listAmps,1));
        for j=1:length(listAmps)
            for i=1:length(listAmps)
                Diff(i,j)=abs(listAmps(i)-listAmps(j));
            end
        end
        
        
        Diff=Diff/max(max(params.arrayInfo.listAmpsAll));
        
        
        params = MakeStimKernels(params,useBrownianStimElec,Diff);
        
        
        
    end
end

if(useKernels)
    if(params.bundle.useBundleAlg)
        [spikes, Log, params]=SpikeSortingAllCases(params,TracesAll,xPostBundle);
    else
        [spikes, Log, params]=SpikeSortingAllCases(params,TracesAll);
    end
else
    [spikes Log params]=SpikeSortingNoKer(params,TracesAll,meanTraces);
end

if(useStimElectrodeAfterBundle+useStimElectrodeBeforeBundle>=1)
    Output.stimInfo.ActiveElectrodes=activeElecs;
    Output.stimInfo.breakpoints=params.patternInfo.breakpoints;
    if(useKernels&&params.global.extraStimElectrode)
        Output.stimInfo.KersStim=params.patternInfo.KersStim;
        Output.stimInfo.xStim=params.patternInfo.xStim;
    end
else
    Output.stimInfo.ActiveElectrodes=activeElecs(ind);
end


Output.neuronInfo=params.neuronInfo;
Output.neuronInfo.neuronIds=neuronIds;
Output.neuronInfo.spikes=spikes;

Output.stimInfo.patternNo=patternNo;
Output.stimInfo.stimElecs=stimElecs;
Output.stimInfo.listAmps=listAmps;
Output.stimInfo.listCurrents=listCurrents;
Output.stimInfo.var0=var0;
Output.stimInfo.channelsConnected=channelsConnected;
if(params.global.saveArt==1)
    Output.stimInfo.Arts=params.patternInfo.Arts;
    Output.stimInfo.firstArt=firstArt;
end
Output.stimInfo.nTrials=params.patternInfo.nTrials;
Output.stimInfo.Kers=params.patternInfo.Kers;

Output.stimInfo.rho=rho;
Output.Log=Log;
try
    Output.path.pathToEi=params.path.pathToEi;
end

Output.path.pathToPreparation=pathToPreparation;
Output.path.pathToAnalysisData=pathToAnalysisData;
Output.stimInfo.sampledTrials=sampledTrials;
Output.bundle=params.bundle;
Output.bundle.onsetC=onsetC;
Output.bundle.onset=onset;
Output.bundle.pvals=pval;
Output.bundle.useStimElectrodeBeforeBundle=useStimElectrodeBeforeBundle;
Output.bundle.useStimElectrodeAfterBundle=useStimElectrodeAfterBundle;
Output.stimInfo.Residuals=Res;
Output.arrayInfo=params.arrayInfo;
Output.stimInfo.polarity=polarity;
Output.arrayInfo.xnew=x(end);
end
