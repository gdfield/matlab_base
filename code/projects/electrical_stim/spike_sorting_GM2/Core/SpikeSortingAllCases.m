function [spikes, Log, params]=SpikeSortingAllCases(params,TracesAll,varargin)
%Gonzalo Mena, 3/2016

Kers=params.patternInfo.Kers;

options=params.global.options;
Q=params.patternInfo.Q;
Qt=params.patternInfo.Qt;
dL=params.patternInfo.dL;
ind=params.patternInfo.indRel;
templates=params.neuronInfo.templates;
Art=params.patternInfo.Art;
saveArt=params.global.saveArt;
if(saveArt)
    Arts{1}=Art;
end
Difs=params.patternInfo.Difs;
stimElec=params.patternInfo.stimElecsRel;
Diags=params.patternInfo.Diags;
useBundleAlg=params.bundle.useBundleAlg;
useStimElectrodeBeforeBundle = params.global.useStimElectrodeBeforeBundle;
useStimElectrodeAfterBundle = params.global.useStimElectrodeAfterBundle;
useBrownian = params.global.useBrownian;
extraStimElectrode     = params.global.extraStimElectrode;
filterStimElectrode    = params.global.filterStimElectrode;
extraNoStimElectrodes  = params.global.extraNoStimElectrodes;
filterNoStimElectrodes = params.global.filterNoStimElectrodes;

elExtra=zeros(size(Art,1),1);

if(useStimElectrodeAfterBundle+useStimElectrodeBeforeBundle>=1)
    if(extraStimElectrode+filterStimElectrode>=1)
        KersSti=params.patternInfo.KersStim;
        
        dLSti=params.patternInfo.dLStim;
        QSti=params.patternInfo.QStim;
        QtSti=params.patternInfo.QtStim;
        xSti=params.patternInfo.xStim;
    end
    breakpoints=params.patternInfo.breakpoints;
    br=intersect(unique(breakpoints{1}(2:end)+1),[1:size(Art,1)]);
    elExtra(br)=stimElec;
else
    
    breakIni=100000;
end

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


nVec=params.bundle.nVec;
updateFreq=params.bundle.updateFreq;
x=params.arrayInfo.x;


els=[];
for n=1:length(templates)
    
    spikes{n}=NaN*zeros(maxCond,size(TracesAll,2));
    [a, b]=sort(max(abs(templates{n}')),'descend');
    ind2=find(a>thresEI);
    if(isempty(ind2))
        ind2=1;
    else
        
    end
    if((length(ind2)==1)&&(b(1)==stimElec))
        ind2=[1 2];
    end
    
    params.neuronInfo.ActiveElectrodes{n}=b(ind2);
    els=union(b(ind2),els);
end

tarray2=setdiff(tarray,0);
params.neuronInfo.ActiveElectrodesAll=els;

if(useStimElectrodeBeforeBundle==0)
    if(isnan(params.bundle.onsBundle))
        elExtra(1:maxCond)=stimElec;
    else
        elExtra(1:params.bundle.onsBundle-1)=stimElec;
    end
end
if(useStimElectrodeAfterBundle==0)
    if(~isnan(params.bundle.onsBundle))
        elExtra(params.bundle.onsBundle:end)=stimElec;
    end
end


indels1=1:Tmax*length(els);
indpat=find(stimElec==els);
indels2=setdiff(1:Tmax*length(els),indpat:length(els):Tmax*length(els));


KnAll=zeros(length(els)*Tmax,1);
indNeurons=[0 kron(1:length(templates),ones(1,length(tarray2)))];
indTimes=[0 kron(ones(1,length(templates)),tarray2)];
indTimesArray=[0 kron(ones(1,length(templates)),1:length(tarray2))];

for n=1:length(templates)
    for t=1:length(tarray2)
        [ActionPotential]=makeActionPotential(n,tarray2(t),templates,Tmax);
        
        Knn(n,t,:,:)=ActionPotential(:,:);
        Aux(:,t)=reshape(ActionPotential(els,:),Tmax*length(els),1)';
        
    end
    KnAll=[KnAll Aux];
end

KnnReshaped=reshape(Knn,size(Knn,1)*size(Knn,2),size(Knn,3),size(Knn,4));

krondiag0=1;
for k=1:2
    krondiag0=kron(krondiag0,dL{k});
end
xold=x;
KersOld=Kers;
Qold=Q;
dLold=dL;
Qtold=Qt;
krondiag0old=krondiag0;

if(useStimElectrodeAfterBundle+useStimElectrodeBeforeBundle>=1)
    indStim=find(1<=breakpoints{1}(2:end));
        
    if(extraStimElectrode+filterStimElectrode>=1)
        indStim=indStim(1);
        indStimOld=indStim;
        dLStim=dLSti{indStim};
        QStim=QSti{indStim};
        QtStim=QtSti{indStim};
        KersStim=KersSti{indStim};
        xStim=xSti(indStim,:);
        
        krondiag0Stim=1;
        krondiag0Stim=kron(krondiag0Stim,dLStim{1});
    end
end
x01=xold(end);

ionset=params.bundle.onsBundle;
if(useBundleAlg==0)
    ionset=1000;
end

try
    nswipes=params.global.nswipes;
catch
    nswipes=1;
end
for nn=1:nswipes
    for i=1:maxCond
        
        if(useStimElectrodeAfterBundle+useStimElectrodeBeforeBundle>=1)
            
                indStim=find(i<=breakpoints{1}(2:end));
                indStim=indStim(1);
                 breakIni=breakpoints{1}(indStim)+1;
           
            
             if(extraStimElectrode+filterStimElectrode>=1)
                
                if(~(indStim==indStimOld))
                    
                    dLStim=dLSti{indStim};
                    QStim=QSti{indStim};
                    QtStim=QtSti{indStim};
                    KersStim=KersSti{indStim};
                    xStim=xSti(indStim,:);
                    
                    krondiag0Stim=1;
                    krondiag0Stim=kron(krondiag0Stim,dLStim{1});
                    
                    indStimOld=indStim;
                end
                
                krondiaginvStim=(exp(xStim(end))*krondiag0Stim*KersStim{2}(i-breakIni+1,i-breakIni+1)+var0).^(-1);
            end
            
        end
        
        
        trialI=nansum(~isnan(squeeze(TracesAll(i,:,1,1))));
        params.patternInfo.nTrials(i)=trialI;
        
        
        els2=setdiff(els,elExtra(i));
        if(elExtra(i)>0)
            indels=indels2;
        else
            indels=indels1;
        end
        
        if(i>=ionset&&ionset>1)
            
            if(i==ionset)
                x=varargin{1};
                if(useBrownian)
                    factp=[0 3 6 7];
                    types=[1 1 8];
                    
                else
                    
                    factp=[0 3 6 9];
                    types=[1 1 1];
                end
                
                for k=1:3
                    [Ker, ~]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),types(k));
                    Kers{k}=Ker;
                    [a, b]=eig(Kers{k});
                    Q{k}=a';
                    Qt{k}=a;
                    dL{k}=diag(b);
                    
                end
                
                
                
                krondiag0=kron(dL{1},dL{2});
               
        
            end
            
            if(mod(i-ionset,updateFreq)==0)
                
                DifAux=Difs;
                DifAux{3}=Difs{3}(i,i);
                DiagsAux=Diags;
                DiagsAux{3}=Diags{3}(i);
                
                if(useBrownian)
                    types=[1 1 8];
                    factp=[0 3 6 7];
                    f1=@(Art,y)logDetKron(Art(i,ind,:),[x(1:7) y(1) log(var0)],DifAux,8,types,DiagsAux,[3 3 1]);
                    g1=@(y)f1(Art,y);
                    x1 = fminunc(g1,[x(end)],params.global.options);
                    x([8])=x1;
                else
                    types=[1 1 1];
                    factp=[0 3 6 9];
                    f1=@(Art,y)logDetKron(Art(i,ind,:),[x(1:9) y(1) log(var0)],DifAux,[10],types,DiagsAux,[3 3 3]);
                    g1=@(y)f1(Art,y);
                    x1 = fminunc(g1,[x(end)],params.global.options);
                    x([10])=x1;
                end
                
                params.arrayInfo.x=x;
                
                
                
                
            end
        end
        
        krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);
        
        
        
        if(i>1)
            if(extraNoStimElectrodes)
                [Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,ArtF(:,ind,:),x,var0);
            else
                [Apred1]=ArtF(i-1,ind,:);
            end
            
            Apred(:,ind,:)=Apred1;
            
            if(size(ArtF,1)<=breakIni||extraStimElectrode==0)
                Apred2=ArtF(end,stimElec,:);
                
            else
                [Apred2]=ExtrapolateArtifactCondStim(KersStim,QStim,QtStim,dLStim,i-breakIni+1,squeeze(ArtF(breakIni:end,stimElec,:)),xStim,var0);
                
            end
            
            Apred(:,stimElec,:)=Apred2;
        else
            Apred=zeros(1,size(Art,2),size(Art,3));
        end
        
        if(saveArt)
            Apreds(i,:,:)=Apred;
        end
        flag=1;
        cont=1;
        
        while(flag==1&&cont<=maxIter)
            
            clear times
            
            if((cont>1)&&(useStimElectrodeAfterBundle+useStimElectrodeBeforeBundle>=1))
                els2=els;
                indels=indels1;
            end
            
            
            AA0=reshape(Apred(:,els2,:),Tmax*length(els2),1);
            
            if(trialI>1)
                TracesResidual=squeeze(TracesAll(i,1:trialI,:,1:Tmax));
            else
                TracesResidual=reshape(squeeze(TracesAll(i,1:trialI,:,1:Tmax)),1,size(TracesAll,3),size(TracesAll,4));
            end
            
            
            
            indTrial=[1:trialI];
            contFindNeurons=1;
            times=zeros(length(templates),trialI);
            
            while(~isempty(indTrial))
                AA=reshape(TracesResidual(:,els2,1:Tmax),trialI,Tmax*length(els2))'-repmat(AA0,1,trialI);
                corrs=-2*AA'*KnAll(indels,:)+repmat(nansum(KnAll(indels,:).^2),trialI,1);
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
            
            if(filterNoStimElectrodes)
                ArtF(i,ind,:)=FilterArtifactLocal(Kers,Art(1:i,ind,:),[x log(var0)],i,ind,Q,Qt,krondiaginv);
            else
                ArtF(i,ind,:)=Art(i,ind,:);
            end
            
            if(~(filterStimElectrode))
                ArtF(i,stimElec,:)=Art(i,stimElec,:);
            else
                
                if(breakIni<i)
                    ArtF(i,stimElec,:)=FilterArtifactLocalStim(KersStim,squeeze(Art(breakIni:i,stimElec,:)),[xStim log(var0)],i-breakIni+1,QStim,QtStim,krondiaginvStim);
                else
                    ArtF(i,stimElec,:)=FilterArtifactLocalStim(KersStim,reshape(squeeze(Art(breakIni:i,stimElec,:)),1,Tmax),[xStim log(var0)],i-breakIni+1,QStim,QtStim,krondiaginvStim);
                end
                
            end
            
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
         
        if((saveArt))
                ArtFiltered(i,:,:)=Apred;
                ArtNoFiltered(i,:,:)=Art(i,:,:);
            end
        Log.Iter(i)=cont-1;
        
    end
end
if(saveArt)
    Arts{2}=ArtFiltered;
    Arts{3}=Apreds;
    Arts{4}=ArtNoFiltered;
    params.patternInfo.Arts=Arts;
end
