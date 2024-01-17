function [TracesAll, Art, var0, listAmps, listCurrents, stimElecs, onset, onsetC, pval, Res,  sampledTrials, Trace0, polarity, firstArt, channelsConnected]=loadTracesArtSort(pathToAnalysisData,patternNo,Tmax,nTrial,subSampleRate,arrayObj,varargin)
%load data from a pattern in a given folder
%also loads stimulus data and construct a firt artifact estimate (means
%accross trials). nTrial is the maximum number of trials
%It authomatically sorts traces increasingly with amplitude
%Gonzalo Mena,3/2016
% Optional inputs. 
%     params


%Set default parameters
Res =NaN;
findBundle=0;
params = []; 
sampledTrials=[];
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
        case 'params'
            params=varargin{j*2}; 
        case 'sampledtrials'
            sampledTrials=varargin{j*2};
        case 'nsubamplitudes'
            nsubAmplitudes=varargin{j*2};    
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end


if(~isempty(params))
    Tmax=params.global.Tmax;
    findBundle=params.bundle.findBundle;
    findBundleTimes=params.bundle.findBundleTimes;
    nNeighborsExclude=params.bundle.nNeighborsExclude;
    detectThreshold= params.bundle.detectThreshold;
end

% Set number of channels to get data from.
activeElectrodes = arrayObj.getElectrodes; 
numActiveElecs = length(activeElectrodes); 
movieNos  = findMovieNos(pathToAnalysisData,patternNo);
[~, channelsWithStim, ~, ~, ~, ~] = getStimAmps(pathToAnalysisData, patternNo, movieNos(1),'numElectrodes',arrayObj.getNumElecs);

% Allocate memory
TracesAll = NaN(length(movieNos),nTrial,numActiveElecs,Tmax);
nTrials=zeros(1,length(movieNos));
listStimElecs=zeros(length(movieNos),length(channelsWithStim)); 


%take entire trace for first amplitude of stimulation
    dataTraces = NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo, movieNos(1), 99999);  
      Trace0=squeeze(dataTraces(unidrnd(size(dataTraces,1)),:,:))-squeeze(mean(dataTraces(:,:,:),1));

for m=1:length(movieNos);
    
    dataTraces = NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo, movieNos(m), 99999);  
    [amps, channelsWithStim, stimAmpVectors, channelsConnected, elecCurrentStep, currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m),'numElectrodes',arrayObj.getNumElecs);
    listStimElecs(m,:)         = channelsWithStim;
    Amp=abs(amps(:,:))';
   polarity=sign(amps(:,:));
    if(m==1)
        nTrials(m)=size(dataTraces,1);
        listAmps(1,:)=Amp;
        listCurrents(1,:)=currentRangesUsed(1,:);
        TracesAll(1,1:size(dataTraces,1),1:numActiveElecs,:)=dataTraces(:,activeElectrodes,1:Tmax);        
        continue
    end
    
   findA=[];
       for k=1:size(listAmps,2)
           if(isempty(find(Amp(:,k)==listAmps(:,k))))
               findA= [];
               break
               
           else
                
               findA=[findA find(Amp(:,k)==listAmps(:,k))'];
           end
       end
       findAmp=[];
       findAU=unique(findA);
       for k=1:length(findAU)
           a=find(findAU(k)==findA);
           if(length(a)==size(listAmps,2))
               findAmp=findAU(k);
           end
       end
    
    if(isempty(findAmp))
        aux=[listAmps(:,1);Amp(:,1)];
        [~, b]=sort([listAmps(:,1);Amp(:,1)],'ascend');
             
        index=find(b==length(aux));
        if(index<length(aux))
            
            TracesAllOld=TracesAll(1:length(listAmps),:,:,:);        
            TracesAll(index,1:size(dataTraces,1),1:numActiveElecs,:)=dataTraces(:,activeElectrodes,1:Tmax);
            
            TracesAll(setdiff([1:length(aux)],index),:,:,:)=TracesAllOld;
            nTrialsOld=nTrials(1:size(listAmps,1));
            nTrials(index)=size(dataTraces,1);
            nTrials(setdiff([1:length(aux)],index))=nTrialsOld;
            
            listAmps=sort([listAmps;Amp],'ascend');
            listCurrentsOld=listCurrents;
            listCurrents(index,:)=currentRangesUsed(1,:);
            listCurrents(setdiff([1:length(aux)],index),:)=listCurrentsOld;
        else
            
            TracesAll(index,1:size(dataTraces,1),1:numActiveElecs,:)=dataTraces(:,activeElectrodes,1:Tmax);
            nTrialsOld=nTrials(1:size(listAmps,1));
            nTrials(index)=size(dataTraces,1);
            nTrials(setdiff([1:length(aux)],index))=nTrialsOld;
            
            listAmps=sort([listAmps;Amp],'ascend');
            listCurrentsOld=listCurrents;
            listCurrents(index,:)=currentRangesUsed(1,:);
            listCurrents(setdiff([1:length(aux)],index),:)=listCurrentsOld;
        end
        
    else
        TracesAll(findAmp,nTrials(findAmp)+1:nTrials(findAmp)+size(dataTraces,1),1:numActiveElecs,:)=dataTraces(:,activeElectrodes,1:Tmax);
        nTrials(findAmp)=nTrials(findAmp)+size(dataTraces,1);
    end
    
end


% Allocate memory

if(isempty(sampledTrials))
    sampledTrials = NaN*zeros(length(listAmps),max(nTrials));
    nTrialsSub = floor(nTrials*subSampleRate);
    indsampleTrials=1;
else
    indsampleTrials=0;
    nTrialsSub=nansum(sampledTrials);
end

stimElecs = unique(listStimElecs);

for e=1:length(stimElecs)
    stimElecsRel(e)=find(activeElectrodes==stimElecs(e));
end

TracesAllOld = TracesAll(1:length(listAmps),1:max(nTrials),:,:);
TracesAll = NaN(length(listAmps),max(nTrialsSub),numActiveElecs,Tmax);


varm = nan(5,1); 
Art = nan(length(listAmps),numActiveElecs,Tmax); 

for m = 1:length(listAmps)
    if(indsampleTrials==1)
    sample=sort(randsample(nTrials(m),nTrialsSub(m)));
    sampledTrials(m,sample)=1;
    else
        sample=find(sampledTrials(m,:)==1);
        nTrialsSub(m)=length(sample);
    end
    if(m==1)
       firstArt=nanmean(TracesAllOld(1,sample,:,:),2);
    end
    
     TracesAll(m,1:nTrialsSub(m),:,:)=TracesAllOld(m,sample,:,:)-repmat(firstArt(:,:,:,:),1,nTrialsSub(m),1);      
    
   if(m<=5)
        a=TracesAll(m,:,setdiff(1:numActiveElecs,stimElecsRel),:);      
        varm(m)=nanvar(a(:));
    end
    Art(m,:,:)=squeeze(nanmean(TracesAll(m,:,:,:),2));   
end

var0=nanmean(varm(1:5));


if(~isempty(nsubAmplitudes))
    listAmps=listAmps(1:nsubAmplitudes:end);
    Art=Art(1:nsubAmplitudes:end,:,:);
    TracesAll=TracesAll(1:nsubAmplitudes:end,:,:,:);
    listCurrents=listCurrents(1:nsubAmplitudes:end);
    if(~isempty(sampledTrials))
        sampledTrials=sampledTrials(1:nsubAmplitudes:end,:);
    end
    firstArt=firstArt(1:nsubAmplitudes:end,:,:,:);
    Trace0=Trace0(1:nsubAmplitudes:end,:,:,:);
end
    
if(findBundle)
    
    [Res]=ResidualsElectrodeSimple(Art,stimElecsRel,findBundleTimes); 
    elsExclude=[];
    for e=1:length(stimElecs);
    elsExclude =[elsExclude arrayObj.getNeighbors(stimElecs(e),nNeighborsExclude)];
    end
    elsExclude=unique(elsExclude);
    nConds=size(Art,1);
    pval = nan(nConds,1); 
    channels = setdiff(activeElectrodes,elsExclude); 
    res_idx = zeros(length(channels),1); 
    for c = 1:length(channels)
        res_idx(c) = find(channels(c)==activeElectrodes);
    end
    
    for k=2:nConds        
        [~, p]=kstest((log(Res(k,res_idx))-nanmean(log(Res(k,res_idx))))./nanstd((log(Res(k,res_idx)))));        
        pval(k)=log(p);
    end
    
    pval(1)=NaN;
    
    aux=find(pval>detectThreshold);
    comp=lastConnectedComponent(aux);
    
    if(isempty(aux));
        onset=listAmps(2);
        onsetC=2;
        return
    else
        if(length(comp)==1&&comp(end)==nConds)
            if(aux(1)==2)
                onset=NaN;
                onsetC=NaN;
                return
            else
                condi=2;
            end
        elseif(length(comp)==1&&comp(end)<nConds);
            condi=comp(end)+1;
        else
            if(comp(end)==nConds)
                condi=comp(end-1)+1;
            else
                condi=comp(end)+1;
            end
        end
        
        onset=listAmps(condi);
        onsetC=condi;
    end
else
    onset=NaN;
    onsetC=NaN;
    pval=NaN;
end


