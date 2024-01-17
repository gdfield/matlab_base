function [Output simParams]=DoSimulateTwoStim(TracesAll,stimElecs,Art,params,listAmps,templates)
%Gonzalo Mena, 03/2016
params.global.Tmax=40;
params.global.tarray=[0 [7:30]];
params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on');

positions=params.global.positions;
x=params.arrayInfo.x;
useBundleAlg = 0;
useStimElectrodeBeforeBundle = 0;
useStimElectrodeAfterBundle = 0;

Tmax=params.global.Tmax;

 

TracesAll2=TracesAll;

%% parameters
pneu =1; %amount of spiking neurons
pspi =[1]; %probability of neurons spiking
xProj = listAmps;
neuronIds=[];

contneu=0;
for n=1:length(templates)
    if(unifrnd(0,1)<pneu)
        contneu=contneu+1;
        neuronIds=[neuronIds n];
        templatesLocal{contneu}=templates{n};
        [a b]=sort(max(abs(templatesLocal{contneu}')),'descend');
       prob=pspi
        if(unifrnd(0,1)<prob)
            
            thresH=unifrnd(0,max(listAmps));
            erfParams(1)=unifrnd(0,4);
            erfParams(2)=-erfParams(1)*thresH;
        else
            thresH=2*max(listAmps);
            erfParams(1)=unifrnd(0,4);
            erfParams(2)=-erfParams(1)*thresH;
        end
        projection = 0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));

        projections(contneu,:)=projection;
            for j=1:length(listAmps)
                nTrials(j)=length(find(~isnan(TracesAll(j,:,1,1))));
            end
            spikesTrue{contneu}=NaN*zeros(length(listAmps),max(nTrials));
            
            for j=1:length(listAmps)
                %spikesTrue{contneu}(j,1:nTrials(j))=(unifrnd(0,1,1,nTrials(j))<projection(j)).*(unidrnd(30,1,nTrials(j))+6);
                
                t0=unidrnd(24)+6;
                spikesTrue{contneu}(j,1:nTrials(j))=(unifrnd(0,1,1,nTrials(j))<projection(j)).*max(7,min(floor(normrnd(t0,25*(1-projection(j)),1,nTrials(j))),30));
                
                for i=1:nTrials(j)
                 [ActionPotential]=makeActionPotential(contneu,spikesTrue{contneu}(j,i),templatesLocal,Tmax);
                 TracesAll2(j,i,:,:)=TracesAll2(j,i,:,:)+reshape(ActionPotential,1,1,512,Tmax);
            end
            end
    end
end
    

params.patternInfo.stimElecs= stimElecs;
params.patternInfo.Art=Art;
params.global.Tmax=40;


%templates=makeTemplatesFromEiShift(pathToEi, neuronIds,[1:512]);
 params.neuronInfo.templates = templatesLocal;

 
 

[spikes Log params]=spikeSortingTwoElectrodesNoStim(params,TracesAll2);
   




Output.neuronInfo.spikes=spikes;
%Output.stimInfo.Art=params.patternInfo.Art;



simParams.spikesTrue=spikesTrue;
simParams.projections=projections;
simParams.neuronIds=neuronIds;
simParams.pneu = pneu; %amount of spiking neurons
simParams.pspi=pspi;
end