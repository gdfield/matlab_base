function [Measures Activations thresH thresAlg  nTrials timesFN params onsets timesAll projections]=CompareElecRespSpike(pathToAnalysisData,neuronId,patternNo,Output,off,thresHolds,varargin)

listAmps=Output.stimInfo.listAmps;
nneu=find(Output.neuronInfo.neuronIds==neuronId);
spikes=Output.neuronInfo.spikes{nneu};
name=['elecResp_n' num2str(neuronId) '_p' num2str(patternNo) '.mat'];
load([pathToAnalysisData name]);
sp=NaN*zeros(size(spikes));
elecStimAmps=abs(elecResp.stimInfo.stimAmps);
nTrialsMax=size(spikes,2);

for i=1:length(listAmps)
    ind=find(listAmps(i)==elecStimAmps);
    
    spaux=[];
    for j=1:length(ind)
        spaux=[spaux elecResp.analysis.latencies{ind(j)}'];
    end
    sp(i,1:length(spaux))=spaux;
   
end


if(nargin==7)
    
    sampledTrials=varargin{1};
else
    sampledTrials=ones(size(spikes,1),size(spikes,2));
end


spaux2=NaN*zeros(size(sp,1),max(nansum(sampledTrials')));

for i=1:size(sampledTrials,1)
spaux2(i,1:nansum(sampledTrials(i,:)==1))=sp(i,sampledTrials(i,:)==1);
end

sp=spaux2;
    
for i=1:size(spikes)
     nTrials(i)=length(find(~isnan(sp(i,:))));
end

Jeff0=length(find(nTrials>0));

sp0=sp;
spikes0=spikes;

for p=1:length(thresHolds)
    [onset onsetC]=findBundleFrompValStandAlone(Output.bundle.pvals,listAmps,thresHolds(p));
    
    onsets(p,:)=[onset onsetC];
    
    
    if(~isnan(onset))
        
        Jeff=min(Jeff0,onsetC-1);
    else
        Jeff=Jeff0;
    end
    
    
    
    spikes=spikes0(1:Jeff,:);
    sp=sp0(1:Jeff,:);
    
    
    spikeCAlg=spikes(:);
    spikeCH=sp(:);
    
    indnonan=find(~isnan(spikeCH));
    spikeCAlg=double(spikeCAlg(indnonan));
    spikeCH=double(spikeCH(indnonan));
    ind1=find(spikeCAlg>0);
    ind2=find(spikeCH>0);
    latCAlg=spikeCAlg(intersect(ind1,ind2));;
    latCH=spikeCH(intersect(ind1,ind2));
    latCH2=spikeCH;
    timesAll{p}=[latCH latCAlg-off];
    
    spikeCH=spikeCH>0;
    spikeCAlg=spikeCAlg>0;
    
    
    Measures(p,:)=[length(intersect(find(spikeCAlg==1),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==0))) length(intersect(find(spikeCAlg==1),find(spikeCH==0)))];
    
    timesFN{p}=latCH2(intersect(find(spikeCAlg==0),find(spikeCH==1)));
    Activations{p}=[nanmean(spikes'>0,1);nanmean(sp'>0,1)];
    [erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff)'; nanmean(spikes'>0,1); nTrials(1:Jeff) ],2, -1);
    projections=projectionComplete;
    thresAlg(p) = -erfParams(2)/erfParams(1);
    params(p,1:2)=erfParams;
    [erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff)'; nanmean(sp(1:Jeff,:)'>0,1); nTrials(1:Jeff) ],2, -1);
    thresH(p)= -erfParams(2)/erfParams(1);
    params(p,3:4)=erfParams;
    projections(3,:)=projectionComplete(2,:);
end


