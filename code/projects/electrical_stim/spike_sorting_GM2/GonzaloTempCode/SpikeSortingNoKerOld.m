function [spikes Log params Apreds]=SpikeSortingNoKer(params,TracesAll,varargin)
%Gonzalo Mena, 3/2016

options=params.global.options;
ind=params.patternInfo.indRel;
templates=params.neuronInfo.templates;
Art=params.patternInfo.Art;
Difs=params.patternInfo.Difs;
Diags=params.patternInfo.Diags;
useBundleAlg=params.bundle.useBundleAlg;
useStimElectrodeBeforeBundle = params.global.useStimElectrodeBeforeBundle;
useStimElectrodeAfterBundle = params.global.useStimElectrodeAfterBundle;
var0=params.patternInfo.var0;
listCurrents=params.patternInfo.listCurrents;
difTrials=1;

thresEI=params.global.thresEI;
Tmax=params.global.Tmax;
tarray=params.global.tarray;
maxIter=params.global.maxIter;
cutBundle=params.bundle.cutBundle;
contMessage=1;
if(cutBundle==1)
    maxCond=params.bundle.onsBundle-1;
else
    maxCond=size(TracesAll,1);
end

stimElecs=params.patternInfo.stimElecsRel;
onsetC=[];
if(params.bundle.findBundle)
    if(~isnan(params.bundle.onsBundle)&&params.bundle.onsBundle>1)
        onsetC=params.bundle.onsBundle-1;
    end
end


breakpoints=findBreakStimElecs(listCurrents);
for e=1:length(stimElecs)
    breakpoints{e}=sort(unique([0 onsetC breakpoints{e}' size(Art,1)]));
end
params.patternInfo.breakpoints=breakpoints;



elExtra=zeros(size(Art,1),length(stimElecs));
for e=1:length(stimElecs)
    br=intersect(unique(breakpoints{e}(2:end)+1),[1:size(Art,1)]);
    elExtra(br,e)=repmat(stimElecs(e),length(br),1);
end

if(useStimElectrodeBeforeBundle==0)
    if(isnan(params.bundle.onsBundle))
        elExtra(1:maxCond,:)=repmat(stimElecs,maxCond,1);
    else
        elExtra(1:params.bundle.onsBundle-1,:)=repmat(stimElecs,params.bundle.onsBundle-1,1);
    end
    
    
end
if(useStimElectrodeAfterBundle==0)
    if(~isnan(params.bundle.onsBundle))
        elExtra(params.bundle.onsBundle:end,:)=repmat(stimElecs,size(Art,1)-params.bundle.onsBundle+1,1);
    end
end



els=[];

for n=1:length(templates)
    
    spikes{n}=NaN*zeros(maxCond,size(TracesAll,2));
    [a b]=sort(max(abs(templates{n}')),'descend');
    elsAux=[];
    for e=1:length(stimElecs)
        indStim(e)=find(stimElecs(e)==b);
        if(a(indStim(e))>thresEI)
            elsAux=[elsAux stimElecs(e)];
        end
    end
    
    ind2=find(a>thresEI);
    inddif=setdiff(ind2,indStim);
    if(isempty(setdiff(ind2,indStim)))
        
        inddif=setdiff([1:length(a)],indStim);
        inddif=inddif(1);
        
    end
    
    
    params.neuronInfo.ActiveElectrodes{n}=union(elsAux,b(inddif));
    els=union(union(elsAux,b(inddif)),els);
end

tarray2=setdiff(tarray,0);
params.neuronInfo.ActiveElectrodesAll=els;



for i=1:size(elExtra,1);
    indels{i}=1:Tmax*length(els);
    indpat=[];
    for e=1:size(elExtra(i,:),2)
        
        indpat=union(indpat,find(elExtra(i,e)==els));
    end
    
    for j=1:length(indpat)
        indels{i}=setdiff(indels{i},indpat(j):length(els):Tmax*length(els));
    end
end







KnAll=zeros(length(els)*Tmax,1);
indNeurons=[0 kron([1:length(templates)],ones(1,length(tarray2)))];
indTimes=[0 kron(ones(1,length(templates)),tarray2)];
indTimesArray=[0 kron(ones(1,length(templates)),[1:length(tarray2)])];

for n=1:length(templates)
    for t=1:length(tarray2)
        [ActionPotential]=makeActionPotential(n,tarray2(t),templates,Tmax);
        
        Knn(n,t,:,:)=ActionPotential(:,:);
        Aux(:,t)=reshape(ActionPotential(els,:),Tmax*length(els),1)';
        
    end
    KnAll=[KnAll Aux];
end

KnnReshaped=reshape(Knn,size(Knn,1)*size(Knn,2),size(Knn,3),size(Knn,4));


flag=1;


i=1;



els2=setdiff(els,elExtra(i,:));

trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
params.patternInfo.nTrials(i)=trialI;
cont=1;
while(flag==1&&cont<=maxIter)
    
    clear times
    
    ArtF(:,ind,:)=Art(1,ind,1:Tmax);
    
    
    ArtF(:,stimElecs,:)=Art(1,stimElecs,1:Tmax);
    
    AA0=reshape(ArtF(i,els2,:),Tmax*length(els2),1);
    
    
    if(trialI>1)
        TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
    else
        TracesResidual=reshape(squeeze(TracesAll(i,1:trialI,:,1:Tmax)),1,size(TracesAll,3),size(TracesAll,4));
    end
    
    indTrial=[1:trialI];
    contFindNeurons=1;
    times=zeros(length(templates),trialI);
    while(length(indTrial)>0)
        AA=reshape(TracesResidual(indTrial,els2,1:Tmax),length(indTrial),Tmax*length(els2))'-repmat(AA0,1,length(indTrial));
        
        
        corrs=-2*AA'*KnAll(indels{i},:)+repmat(nansum(KnAll(indels{i},:).^2),length(indTrial),1);
        if(contFindNeurons>1)
            for trial=1:length(indTrial)
                corrs(trial,find(indNeurons==indNeurons(tmax(indTrialRel(trial)))))=NaN;
                
            end
        end
        [mins tmax]=nanmin(corrs');
        
        indTrialRel=find(tmax>1);
        indTrial=indTrial(indTrialRel);
        idx = sub2ind(size(times), indNeurons(tmax(indTrialRel)),indTrialRel);
        times(idx)=indTimes(tmax(indTrialRel));
        
        idxSubtract = sub2ind(size(squeeze(Knn(:,:,1,1))), indNeurons(tmax(indTrialRel)),indTimesArray(tmax(indTrialRel)));
        
        TracesResidual(indTrial,:,:)=TracesResidual(indTrial,:,:)-KnnReshaped(idxSubtract,:,:);
        contFindNeurons=contFindNeurons+1;
    end
    
    
    Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
    
    ArtF(1,ind,:)=Art(1,ind,1:Tmax);
    
    ArtF(1,stimElecs,:)=Art(1,stimElecs,:);
    
    flag2=ones(length(templates),1);
    for n=1:length(templates)
        
        if(nansum(times(n,:)==spikes{n}(1,1:trialI))==trialI)
            flag2(n)=0;
        end
    end
    flag=max(flag2);
    for n=1:length(templates)
        spikes{n}(1,1:trialI)=times(n,:);
        
    end
    cont=cont+1;
end
Log.Iter(1)=cont;





ionset=params.bundle.onsBundle;
if(useBundleAlg==0)
    ionset=1000;
end

for i=2:maxCond
    
    
    els2=setdiff(els,elExtra(i,:));
    
    
    
    trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
    params.patternInfo.nTrials(i)=trialI;
    
    if(isempty(varargin{1}))
        Apred=ArtF(i-1,:,:);
    else
        Apred=zeros(1,size(Art,2),size(Art,3));
    end
    
    Apreds(i-1,:,:)=Apred;
    
    flag=1;
    cont=1;
    while(flag==1&&cont<=maxIter)
        
        
        clear times
        
        AA0=reshape(Apred(:,els2,:),Tmax*length(els2),1);
        
        
        if(trialI>1)
            TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
        else
            TracesResidual=reshape(squeeze(TracesAll(i,1:trialI,:,1:Tmax)),1,size(TracesAll,3),size(TracesAll,4));
        end
        indTrial=[1:trialI];
        contFindNeurons=1;
        times=zeros(length(templates),trialI);
        %indTrialSpike=indTrial;
        while(~isempty(indTrial))
            AA=reshape(TracesResidual(:,els2,1:Tmax),trialI,Tmax*length(els2))'-repmat(AA0,1,trialI);
            
            corrs=-2*AA'*KnAll(indels{i},:)+repmat(nansum(KnAll(indels{i},:).^2),trialI,1);
            if(contFindNeurons>1)
                for trial=1:length(indTrial)
                    corrs(indTrial(trial),find(indNeurons==indNeurons(tmax(indTrial(trial)))))=NaN;
                    
                end
            end
            [mins, tmax]=nanmin(corrs');
            
            indTrial=find(tmax>1);
            
            idx = sub2ind(size(times), indNeurons(tmax(indTrial)),indTrial);
            times(idx)=indTimes(tmax(indTrial));
            
            idxSubtract = sub2ind(size(squeeze(Knn(:,:,1,1))), indNeurons(tmax(indTrial)),indTimesArray(tmax(indTrial)));
            
            TracesResidual(indTrial,:,:)=TracesResidual(indTrial,:,:)-KnnReshaped(idxSubtract,:,:);
            contFindNeurons=contFindNeurons+1;
        end
        
        
        
        
        Art(i,:,:)=squeeze(nanmean(TracesResidual,1));
        ArtF(i,:,:)=Art(i,:,:);
        Apred=ArtF(i,:,:);
        
        
        flag2=ones(length(templates),1);
        for n=1:length(templates)
            
            if(nansum(times(n,:)==spikes{n}(i,1:trialI))>=trialI-difTrials)
                flag2(n)=0;
            end
        end
        flag=max(flag2);
        for n=1:length(templates)
            spikes{n}(i,1:trialI)=times(n,:);
            
        end
        cont=cont+1;
        if(cont==maxIter)
            Log.Message{contMessage}=['Maximum number of iterations exceeded at conditon ' num2str(i) ];
            contMessage=contMessage+1;
        end
    end
    
    Log.Iter(i)=cont-1;
    
end

params.patternInfo.Art=Art(1:maxCond,:,:);
