function [Measures Activations thresH thresAlg  timesFN params onsets timesAll]=CompareElecRespSpikeTwoFolders(pathToAnalysisData1,pathToAnalysisData2,pathCollapse,patternNo,varargin)


dirs1=dir(pathToAnalysisData1);
dirs2=dir(pathToAnalysisData2);
str=['p' num2str(patternNo) '.mat'];

for i=1:length(dirs1)
    if(~isempty(strfind(dirs1(i).name,str)))
    
    load([pathToAnalysisData1 dirs1(i).name])
    neuron=elecResp.cells.main;
    aux=elecResp.analysis.latencies;
    
    [listAmpsNew listStimElecs TracesAll4 Art4 channelsConnected]=loadAmps(pathCollapse,patternNo);
    [~,~,indexCollapse]=collapseAmplitudesEqual(listAmpsNew,listAmpsNew);
    for k=1:length(indexCollapse)
        sind(k)=length(indexCollapse{k});
    end
    
    spikesGT=NaN*zeros(length(indexCollapse),length(aux{1})*max(sind));
    for k=1:length(indexCollapse)
        for h=1:length(indexCollapse{k})
            spikesGT(k,(h-1)*length(aux{1})+1:(h-1)*length(aux{1})+length(aux{indexCollapse{k}(h)}))=aux{indexCollapse{k}(h)}';
        end
    end
    
    for j=1:length(dirs2)
        if(~isempty(strfind(dirs2(j).name,str)))
            load([pathToAnalysisData2 dirs2(j).name])
            indneuron=find(elecRespAuto.neuronInfo.neuronIds==neuron);
            spikesAlg=elecRespAuto.neuronInfo.spikes{indneuron};
            onsetC=elecRespAuto.bundle.onsetC;
            
            break
        end
        
    end
    end
end

nTrials=elecRespAuto.stimInfo.nTrials;
listAmps=elecRespAuto.stimInfo.listAmps;

if(size(listAmps,2)>1)
    listAmps=sqrt(nansum(listAmps'.^2))';
end
if(nargin>=5)
    
    thresHolds=varargin{1};
else
    thresHolds=-1000000000;
end

for p=1:length(thresHolds)
    [onset onsetC]=findBundleFrompValStandAlone(elecRespAuto.bundle.pvals,listAmps(:,1),thresHolds(p));
    onsets(p,:)=[onset onsetC];
    
    
    if(~isnan(onset))
        
        Jeff=onsetC-1;
    else
        Jeff=size(spikesAlg,1);
    end
    
    
    
    spikesAlg=spikesAlg(1:Jeff,:);
    spikesGT=spikesGT(1:Jeff,:);
    
    
    spikeCAlg=spikesAlg(:);
    spikeCH=spikesGT(:);
    
    indnonan=find(~isnan(spikeCH));
    spikeCAlg=double(spikeCAlg(indnonan));
    spikeCH=double(spikeCH(indnonan));
    latCAlg=spikeCAlg;
    spikeCAlg=spikeCAlg>0;
    
    latCH=spikeCH;
    
    spikeCH=spikeCH>0;
    off=varargin{2};
    timesAll{p}=[latCH latCAlg-off];
    
    Measures(p,:)=[length(intersect(find(spikeCAlg==1),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==0))) length(intersect(find(spikeCAlg==1),find(spikeCH==0)))];
    
    timesFN{p}=latCH(intersect(find(spikeCAlg==0),find(spikeCH==1)));
    Activations{p}=[nanmean(spikesAlg'>0);nanmean(spikesGT'>0)];
    [erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff)'; nanmean(spikesAlg'>0); nTrials(1:Jeff) ],2, -1);
    thresAlg(p) = -erfParams(2)/erfParams(1);
    params(p,1:2)=erfParams;
    [erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff)'; nanmean(spikesGT(1:Jeff,:)'>0); nTrials(1:Jeff) ],2, -1);
    thresH(p)= -erfParams(2)/erfParams(1);
    params(p,3:4)=erfParams;
end

