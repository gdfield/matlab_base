 pat2=1114;
 pat11=390;
 pat473=424;
 [listAmps3 listStimElecs3 TracesAll3 Art3]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data007',1114);
 
  
 %%Two opposite patterns on electrode 468
 
 pathToPreparation='/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/'
 folders={'data007'};
 
 params=InitializeArray(pathToPreparation,1,folders,patternNos(1));
  patternNos(2)=3789;
  params=InitializeArray(pathToPreparation,1,folders,patternNos(2));
 
  %%Two single electrode stimulation
   pathToPreparation='/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/'
 folders={'data011'};
 
  patternNos(1)=27;
  [listAmps listStimElecs TracesAll1 Art1 pat ch]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data011',27);
  params=InitializeFromArt(1.5*Art1(:,:,1:40),patternNos(1),abs(listAmps),40,20,[],[],[1:10]);
  
  paramsAll(1)=params;
  patternNos(2)=14;
  [listAmps2 listStimElecs2 TracesAll2 Art2 pat ch]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data011',14);
   params=InitializeFromArt(-1.5*Art2(:,:,1:40),patternNos(1),abs(listAmps),40,20,[],[],[1:10]);
  
  
  %% two electrodes, from actual pattern. Fit individual 
  i=450
  patternNos(1)=list(i,11);
  patternNos(2)=list(i,12);
 electrodes=list(i,8:9);
  [listAmps listStimElecs TracesAll Art pat ch]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data007',list(i,10));
 [listAmps1 listStimElecs1 TracesAll1 Art1  ch]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data007', patternNos(1));
 [listAmps2 listStimElecs2 TracesAll2 Art2 ch]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data007', patternNos(2));

 
  x=[0.1890    2.3229    0.7885   -1.8167    0.4993  -30.4399   -1.1175 -100.0000 -100.0000   30.8661];
 
 
  permutation(1:7)=[1:7];
permutation(8:9)=[9:10];
permutation(10)=8;

  params=InitializeFromArt(Art1(:,:,1:40),electrodes(1),abs(listAmps1),40,20,[],[],permutation);
  paramsAll(1)=params;
  
   permutation(1:9)=[2:10];
  permutation(10)=[1];

  params=InitializeFromArt(Art2(:,:,1:40),electrodes(2),abs(listAmps1),40,20,paramsAll(1).arrayInfo.x(1:7),[1:7],permutation);
  paramsAll(2)=params;
  
 
  

KerAux1=NaN*zeros(512,512);
%KerAux1(setdiff([1:512],electrodes(1)),setdiff([1:512],electrodes(1)))=paramsAll(1).patternInfo.Kers{2}*exp(paramsAll(1).arrayInfo.x(end));
KerAux1(setdiff([1:512],electrodes(1)),setdiff([1:512],electrodes(1)))=paramsAll(1).patternInfo.Kers{2};

KerAux1(electrodes(2),:)=NaN;
KerAux1(:,electrodes(2))=NaN;

KerAux2=NaN*zeros(512,512);
KerAux2(setdiff([1:512],electrodes(2)),setdiff([1:512],electrodes(2)))=paramsAll(2).patternInfo.Kers{2}*exp(paramsAll(2).arrayInfo.x(end)-paramsAll(1).arrayInfo.x(end));
%KerAux2(setdiff([1:512],electrodes(2)),setdiff([1:512],electrodes(2)))=paramsAll(2).patternInfo.Kers{2};

KerAux2(electrodes(1),:)=NaN;
KerAux2(:,electrodes(1))=NaN;

Kers{1}=paramsAll(1).patternInfo.Kers{1};
Kers{3}=paramsAll(1).patternInfo.Kers{3};
x=paramsAll(1).arrayInfo.x;

Kers{2}=KerAux1+KerAux2;
Kers{2}=Kers{2}(setdiff([1:512],electrodes),setdiff([1:512],electrodes));


for k=1:3
[a b]=eig(Kers{k});
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
end

params2=paramsAll(1);
params2.patternInfo.Q=Q;
params2.patternInfo.Qt=Qt;
params2.patternInfo.Kers=Kers;
params2.patternInfo.dL=dL;
params2.patternInfo.ind=setdiff([1:512],electrodes);
params2.patternInfo.var0=10;


[Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,17,Art(1:16,setdiff([1:512],electrodes),1:40),x,params2.patternInfo.var0);   
ArtExtra=NaN*zeros(1,512,40);
ArtExtra(1,setdiff([1:512],electrodes),:)=Apred1;
PlotArtifact(Art,17,electrodes(1),electrodes,18,14,1)
PlotArtifact(ArtExtra,1,electrodes(1),electrodes,18,14,2)
PlotArtifact(ArtExtra-Art(17,:,1:40),1,electrodes(1),electrodes,18,14,3)
PlotArtifact(Art,16,electrodes(1),electrodes,18,14,4)
PlotArtifact(Art(16,:,1:40)-Art(17,:,1:40),1,electrodes(1),electrodes,18,14,5)

for j=1
factp=[0 3 6 9];
for j=1:2
Diags=paramsAll(j).patternInfo.Diags;
Difs=paramsAll(j).patternInfo.Difs;
x=paramsAll(j).arrayInfo.x;
rho=Diags{2};
for k=1:3
[Ker KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),[1 1 1]);
[a b]=eig(Ker);
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
Kers{k}=Ker;
Kernels{j,k}=Kers{k};
end
end



  templates{1}=templatesAll{30}
for i=1:100
[Output simParams]=DoSimulateTwoStim(TracesAll(:,:,:,1:40),electrodes,5*Art(:,:,1:40),params2,abs(listAmps1),templates);
err(i)= nansum(nansum(abs(Output.neuronInfo.spikes{1}-simParams.spikesTrue{1})>0))/prod(size(Output.neuronInfo.spikes{1}))
end


for i=1:93
[a b]=sort(max(abs(templatesAll{i}')),'descend');
elgood(i)=b(1);
end