preparations={'2012-09-24-3','2014-09-10-0','2014-11-05-3','2014-11-05-8','2014-11-24-2','2015-04-09-2','2015-04-14-0','2015-05-27-0'};



path{1}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/';
dire{1}={'data003','data004','data005','data006'};
path{2}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-09-10-0/';
dire{2}={'data003'};
path{3}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-05-3/';
dire{3}={'data003','data001','data004'};
path{4}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-05-8/';
dire{4}={'data002','data003'};
path{5}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-24-2/';
dire{5}={'data002'};
path{6}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-04-09-2/';
dire{6}={'data002'};
path{7}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-04-14-0/';
dire{7}={'data001'};
path{8}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-05-27-0/';
dire{8}={'data001'};
dire2=dire;
dire2{6}={'data002','data003'};

modalities=[[0 0];[0 1];[1 0];[1 1]];

for j=[1:4]

for g=setdiff([1:7],6)
    try
tic
[params]=InitializeArrayRobust(path{g},arrayObj,'foldernames',dire{g},'AvoidBundleHyper',modalities(j,1),'whiten',modalities(j,2),'typevarestimate',1);
times1(g)=toc;

paramsModalities1{j}(g)=params;
[params]=InitializeArrayRobust(path{g},arrayObj,'foldernames',dire{g},'AvoidBundleHyper',modalities(j,1),'whiten',modalities(j,2),'typevarestimate',2);

paramsModalities2{j}(g)=params;
end
 end
end

%% SPIKE SORTING

    
    
        for j=[1:4]
for g=setdiff([1:7],6)
        
        paramsModalities1{j}(g).bundle.findBundle=1;
        
        paramsModalities1{j}(g).global.sortData=1;
        paramsModalities1{j}(g).global.nTrial=80;
         paramsModalities1{j}(g).global.subSampleRate=1;

paramsModalities1{j}(g).global.tarray=[7:40];

paramsModalities2{j}(g).global.tarray=[7:40];
          tic
         tic
         
         
        paramsModalities2{j}(g).global.saveArt=0;
         
         
        paramsModalities1{j}(g).bundle.useBundleAlg=1;
        paramsModalities1{j}(g).global.useStimElec=1;
        paramsModalities1{j}(g).global.useStimElectrodeBeforeBundle=1;
        paramsModalities1{j}(g).global.useStimElectrodeAfterBundle=0;
        paramsModalities2{j}(g).bundle.findBundle=1;
        paramsModalities1{j}(g).global.saveArt=0;
        
        
        paramsModalities2{j}(g).global.sortData=1;
        paramsModalities2{j}(g).global.nTrial=80;
         paramsModalities2{j}(g).global.subSampleRate=1;
          tic
         tic
        paramsModalities2{j}(g).bundle.useBundleAlg=1;
       paramsModalities2{j}(g).global.useStimElec=1;
        paramsModalities2{j}(g).global.useStimElectrodeBeforeBundle=1;
        paramsModalities2{j}(g).global.useStimElectrodeAfterBundle=0;
        
    end
    
    end


for g=setdiff([2:7],6)
   %for g=6
    %for g=1
    %to avoid undesided 'aliasing' effects.
   
    
    neuronIds=neuronswr{g};
    if(g==6)
        a1=intersect(find(nstimElecwr==1),find(groupwr==6));
        pat=unique(patternwr(indwr{g},1));
    else
        pat=unique(patternwr(indwr{g},1));
    end
   
    
    
    %for p=1:length(pat)
     for p=1:length(pat)   
for j=[1:4]
     try
        patternNo=pat(p);
       
        %if(g==6)
       % [Output]=DoSpikeSortingLargeScale(path{g},patternNo,neuronIds,paramsAllnew(g),templates2{g},arrayObj,'FoldersNames',dire2{g});
  
        % Respnew{g}{p}=Output;
        
        
   [Output]=DoSpikeSortingLargeScale(path{g},patternNo,neuronIds,paramsModalities1{j}(g),templates2{g},arrayObj,'FoldersNames',dire2{g});
   
          RespModalities1{j,g}{p}=Output;
      
 



   [Output]=DoSpikeSortingLargeScale(path{g},patternNo,neuronIds,paramsModalities2{j}(g),templates2{g},arrayObj,'FoldersNames',dire2{g});
   
          RespModalities2{j,g}{p}=Output;
     catch
         [g p]
     end
     
    end
end
end

all until 34

    
for p=1:99
Resp{j,g}{p}=RespModalities1{j,g}{p};
end
%% Load Raw acuracies per preparation
%g==1
%Outputs{120}.neuronInfo.neuronIds=Outputs{121}.neuronInfo.neuronIds
for g=setdiff([1:7],6)
    
    thresHolds=[-200 -12];

    for j=[1:4]  
 Outputs{2*j-1}= RespModalities1{j,g};
    Outputs{2*j}= RespModalities2{j,g};
end

%Outputs{3}=Resp1WhitenBundle{g};
 %   Outputs{4}=Resp2WhitenBundle{g};
%Outputs{5}=Resp1NoWhitenBundle{g};
 %   Outputs{6}=Resp2NoWhitenBundle{g};
%Outputs{7}=Resp1NoWhitenNoBundle{g};
 %   Outputs{8}=Resp2NoWhitenNoBundle{g};
%Outputs{9}=RespOld{g}
    %Outputs{3}=Resp2Whiten{g};
    clear nans
    clear cont
     for k=[1:8]
    
    pat=unique(patternwr(indwr{g},1));
    a1=indwr{g};
    if(g==6)
        neuronIds=Outputs{k}{30}.neuronInfo.neuronIds;
    else
        neuronIds=Outputs{k}{1}.neuronInfo.neuronIds;
    end
    cont=[];
    
    %for i=1:length(a1)
     for i=1:length(a1)
   i
        try
            neuron=neuronIdwr(a1(i));
            nneu=find(neuron==neuronIds);
            patternNo=patternwr(a1(i));
            patindex=find(patternNo==pat);
            pathToAnalysisData=[Outputs{k}{patindex}.path.pathToAnalysisData ];
            
            listAmps=Outputs{k}{patindex}.stimInfo.listAmps;
           
            [Measures Activation thresH thresAlg  nTrial timeFN paramerf onset]=CompareElecRespSpike(pathToAnalysisData,neuron,patternNo,Outputs{k}{patindex},thresHolds);
            
            for p=1:length(thresHolds)
                measures{k}{g,p}(i,:)= Measures(p,:);
                nTrials{k}{g}(i,p)=nansum(Measures(p,:));
                thres{k}{g,p}(i,:)=[  thresH(p) thresAlg(p)  max(listAmps) paramerf(p,:)];
                Activations{k}{g,p}{i}=Activation{p};
                timesFN{k}{g,p,i}=timeFN{p};
                onsets{k}{g}(i,p)=onset(p,1);
                onsetsC{k}{g}(i,p)=onset(p,2);
            end
            %[Measures Activations thresH thresAlg  nTrial timesFN]=CompareElecRespSpike(pathToAnalysisData,neuron,patternNo,Outputs{patindex},thresholds);
        catch [k g i]
        end
    end
    for p=1:length(thresHolds)
        indZero=find(nansum(measures{k}{g,p}')==0);
        measures{k}{g,p}(indZero,:)=NaN;
        nTrials{k}{g}(indZero,p)=NaN;
        thres{k}{g,p}(indZero,:)=NaN;
    end
     end
end

%% Merge 8 preparations
for k=[1:8]
for p=1:length(thresHolds)
    measuresA{k}{p}=[];
    nTrialsA{k}{p}=[];
    thresA{k}{p}=[];
    contsA{k}{p}=[];
    timesFNAllA{k}{p}=[];
    for g=setdiff([1:7],6)
        
        timesFNAll{k}{g,p}=[];
        cont(g,p)=0;
        for i=1:size(measures{k}{g,1},1)
            cont(g,p)= cont(g,p)+length(union(find(timesFN{k}{g,p,i}>40),find(timesFN{k}{g,p,i}<10)));
            conts{k}{g}(i,p)=length(union(find(timesFN{k}{g,p,i}>40),find(timesFN{k}{g,p,i}<10)));
            timesFNAll{k}{g,p}=[timesFNAll{k}{g,p} timesFN{k}{g,p,i}'];
        end
        timesFNAllA{k}{p}=[timesFNAllA{k}{p} timesFNAll{k}{g,p}];
        
        measuresA{k}{p}=[measuresA{k}{p};measures{k}{g,p}];
        nTrialsA{k}{p}=[nTrialsA{k}{p} nTrials{k}{g}(:,p)'];
        thresA{k}{p}=[thresA{k}{p}; thres{k}{g,p}];
        contsA{k}{p}=[contsA{k}{p} conts{k}{g}(:,p)'];
    end
end
end

%% Overal AnLYSIS. fROM HERE TAKE OVERALL MEASURES
%1 True Positive
% 2 False Negative
%3 True Negative
% False Positive
for k=1:3
for p=1:length(thresHolds)
    measuresA2{k}{p}=measuresA{k}{p};
    
    measuresA2{k}{p}(:,1)=measuresA{k}{p}(:,1)+contsA{k}{p}';
    %measuresA2{k}{p}(:,1)=measuresA{k}{p}(:,1);
    
    
    measuresA2{k}{p}(:,2)=measuresA{k}{p}(:,2)-contsA{k}{p}';
    %measuresA2{p}(:,2)=measuresA{p}(:,2);
    AccA(k,p)=nansum(nansum(measuresA2{k}{p}(:,[1 3])'))./nansum(nansum(measuresA2{k}{p}'));
    TPA(k,p)= nansum(nansum(measuresA2{k}{p}(:,[1 2])'));
    TTA(k,p)=nansum(nansum(measuresA2{k}{p}(:,[1 2 3 4])'));;
    TPRA(k,p)=nansum(nansum(measuresA2{k}{p}(:,[1 2])'))./nansum(nansum(measuresA2{k}{p}'));
    SensA(k,p)=nansum(nansum(measuresA2{k}{p}(:,[1])'))./nansum(nansum(measuresA2{k}{p}(:,[1 2])'));
    SpecA(k,p)=nansum(nansum(measuresA2{k}{p}(:,[3])'))./nansum(nansum(measuresA2{k}{p}(:,[3 4])'));
    prob(k,p)=nansum(nansum(measuresA2{k}{p}(:,[1:4])));
end
end

