function [listAmps listStimElecs TracesAll Art ch Art1 Art2 listCurrents]=loadAmps(pathToAnalysisData,patternNo,varargin)
%load data from a pattern in a given folder
%also loads stimulus data and construct a firt artifact estimate (means
%accross trials). nTrial is the maximum number of trials
%It authomatically sorts traces increasingly with amplitude
%Gonzalo Mena,3/2016
Tmax=100;
nTrial=50;
movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);


TracesAll=NaN*zeros(length(movieNos),nTrial,512,Tmax);

nTrials=zeros(1,length(movieNos));
trialmax=1;
ch=[];
for m=1:length(movieNos);
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
        movieNos(m), 99999);
    
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
      listStimElecs(m,:)         = channelsWithStim;
   listAmps(m,:)=amps;
   TracesAll(m,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
        trialmax=max(trialmax,size(dataTraces,1));
        ch=[ch channelsConnected];
        if(~isequal(channelsConnected,patternNo))
            pat=0;
        end
        currents(m,:)=currentRangesUsed;
        
end
ch=unique(channelsConnected);
TracesAll=TracesAll(:,1:trialmax,:,:);

Art=squeeze(nanmean(TracesAll(:,1:end,:,:),2));
Art1=Art(1,:,:);

Art2=Art-repmat(Art1,size(Art,1),1,1);
%for j=1:size(TracesAll,1)
%TracesAll(j,:,:,:)=TracesAll(j,:,:,:)-reshape(repmat(Art1,trialmax,1,1),1,trialmax,size(Art,2),size(Art,3));
%end
if(~isempty(varargin))

   Art=Art2;
   for j=1:size(TracesAll,1)
TracesAll(j,:,:,:)=TracesAll(j,:,:,:)-reshape(repmat(Art1,trialmax,1,1),1,trialmax,size(Art,2),size(Art,3));
   end
end

listCurrents=currents;