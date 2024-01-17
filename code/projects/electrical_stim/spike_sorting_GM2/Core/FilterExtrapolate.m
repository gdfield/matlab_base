
load('/Users/Cybele/Documents/BigData/RespJournal.mat')

%%%% initial params
pathToPreparationInitial='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/';
params=InitializeArray(pathToPreparationInitial,1,{'data006'},474);



pathToPreparation{1}='/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/';
dire{1}={'data008','data009','data011'};
pathToPreparation{2}='/Volumes/DOS/TTX/2015-03-09-0/';
dire{2}={'data018'};
pathToPreparation{3}='/Volumes/DOS/TTX/2014-08-13-0/';
dire{3}={'data011'};
pathToPreparation{4}='/Volumes/DOS/TTX/2014-11-05-3/';
dire{4}={'data008'};

patterns{1}=[[1:319],[321:384]];
patterns{2}=[129:256];
patterns{3}=[385:512];
patterns{4}=[129:256];

arrayObj=Array(501);
for g=1:1
    [paramsM]=InitializeArrayRobust(pathToPreparation{g},arrayObj,'foldernames',dire{g});
    [paramsM2]=InitializeArrayRobust(pathToPreparation{g},arrayObj,'foldernames',dire{g},'useBrownian',0,'useLocalization',[1 0]) 
    params(g,1)=paramsM;
    params(g,2)=paramsM2;
    for k=1:2
    params(g,k).bundle.findBundle=1;
    
    params(g,k).global.sortData=1;
    params(g,k).global.nTrial=80;
    params(g,k).global.subSampleRate=1;
    
    params(g,k).global.tarray=[7:40];
    params(g,k).bundle.useBundleAlg=0;
    
    params(g,k).global.useStimElectrodeBeforeBundle=0;
    params(g,k).global.useStimElectrodeAfterBundle=0;
    params(g,k).bundle.findBundle=0;
    params(g,k).global.saveArt=0;
    params(g,k).global.thresEI=15;
    params(g,k).global.nswipes=4;
    end
end






typep=2;
index=[2 2 2];
paramsM2=parameters{index(1),index(2),index(3)};
x=paramsM2.arrayInfo.x

patternNo=7;
stimElecs=7
%% Load Matrices, blabla
Tmax=55;
activeElecs = arrayObj.getElectrodes;
numActiveElecs = length(activeElecs);
positions=arrayObj.getPositions;


[listAmps listStimElecsn TracesAll Art patn chn]=loadAmps('/Volumes/MAC OS/Research/EJBigData/EJ-2014-11-05-Processed/data011/',stimElecs(1),1);
for i=1:size(Art,1)
Art(i,:,:)=Art(i,:,:)+reshape(normrnd(0,simParamsIn.noise,size(Art,2),size(Art,3)),1,size(Art,2),size(Art,3));
end

listAmps=abs(listAmps);
[~,rho] = cart2pol(positions(:,1)-positions(stimElecs(1),1),positions(:,2)-positions(stimElecs(1),2));



ind=setdiff(activeElecs,find(rho==0));
if(length(stimElecs)>1)
    ind=setdiff(activeElecs,stimElecs);
end


for e=1:length(stimElecs)
    stimElecsRel(e)=find(activeElecs==stimElecs(e));
end

Dif1= zeros(Tmax);
for i=1:Tmax
    for j=1:Tmax
        Dif1(i,j)=abs(i-j);
    end
end
useBrownian=logical(index(1)-1)
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
Difs{3}=Difs{3}/max(max(paramsM2.arrayInfo.listAmpsAll));
Dif1=Dif1/Tmax;
Difs{1}=Dif1;


Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/arrayObj.maxR;
Diags{3}=listAmps(:,1)/max(paramsM2.arrayInfo.listAmpsAll);
if(useBrownian)
    types=[1 1 8];
    factp=[0 3 6 7];
     nvar=[3 3 1];
    type=[1 1 8];
  
else
    types=[1 1 1];
    factp=[0 3 6 9];
    nvar=[3 3 3];
type=[1 1 1];

end

for m=1:36
    if(m<=36)
        a=Art(m,setdiff(1:numActiveElecs,stimElecsRel),:);
        varm(m)=nanvar(a(:));
    end
end
var0=1;

% %% useless
% %usebrownian
% %             f1=@(Art,x)logDetKron(Art(:,:,:),[x(1:6) 0  x(end) log(var0)],Difs,setdiff([1:8],[7]),type,Diags,nvar);
% %             g1=@(x)f1(Art,x);
% %             x1 = fminunc(g1,[2 0 0 2 0 0 10],options);
% %             x=[x1(1:6) 0 x1(1)];
% %        
%   % ind2=setdiff([1:512],316);
%    
%    %     var0=nanmean(varm);
%     %    var0old=var0;
%      %   onsetC=size(Arts,1)+1;
%       %  DifAux=Difs;
% %         DifAux{3}=Dif3(1:onsetC-1,1:onsetC-1)./max(params(1,1).arrayInfo.listAmpsAll);
% %         DiagsAux=Diags;
% %         DiagsAux{3}=listAmps(1:onsetC-1)/max(params(1,1).arrayInfo.listAmpsAll);
% %         if(useBrownian)
% %             f1=@(Art,y)logDetKron(Art(1:onsetC-1,ind,1:Tmax),[x(1:7) y(1) log(var0old)],DifAux,8,types,DiagsAux,[3 3 1]);
% %             g1=@(y)f1(Art,y);
% %             x1 = fminunc(g1,[x(end)],params(1,1).global.options);
% %             x([8])=x1;
% %         else
% %             f1=@(Arts,y)logDetKron(Arts(1:onsetC-1,:,1:Tmax),[y(1:7) -100 -100 10 log(15)],DifAux,setdiff([1:9],[7:8]),types,DiagsAux,[3 3 3]);
% %             g1=@(y)f1(Arts,y);
% %             x1 = fminunc(g1,[-2 0 0 -2 0 0 -2 ],paramsM2.global.options);
% %             %x([10])=x1;
% %         end
% %         
%         
%          f1=@(Arts,y)logDetKron(Arts(1:onsetC-1,:,1:Tmax),[x1(1:7) -100 -100 y(1) log(15)],DifAux,10,types,DiagsAux,[3 3 3]);
%          g1=@(y)f1(Arts,y);
%          x2 = fminunc(g1,[25],paramsM2.global.options);
%          for i=1:40
%              eval(i)=g1(i);
%          end
%         array=[-5:20];
% 
%        for i=1:length(array)
%            i
%            for j=1:length(array)
%                for k=1:length(array)
%                matrix(i,j)=g1([array(i),array(j)]);
%            end
%            end
%        end
%        
%      A=matrix;
%     [M,I] = min(A(:));
%     [I_row, I_col] = ind2sub(size(A),I)
%        
%        
%         
% i
% 
% %x(7)=-4
% 
%    %x=[10 -100 -100 10 -100 -100 x1(1) -100 -100 x1(2)]       
% x=[x1(1:7) -100 -100 x2];
%           var0=15;
%           
;

%% useful
for k=1:3
    [Ker, KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),types(k));
%   if(k<=2)
%     Ker=eye(size(Difs{k},1));
%    end
     Kers{k}=Ker;
    [a, b]=eig(Kers{k});
    Q{k}=a';
    Qt{k}=a;
    dL{k}=diag(b);
    
end


% f1=@(Art,x)logDetKron(Art(:,ind,1:Tmax),[x(1) -100 -100 x(2) -100 -100  0 x(end-1) x(end)],Difs,setdiff([1:9],[2 3 5 6 7]),types,Diags,nvar);
 %           g1=@(x)f1(Art,x);
  %          x1 = fminunc(g1,[0 0 10 log(var0)],options);
   %         x=[x1(1) -100 -100 x1(2) -100 -100  0 x1(end-1)];
       
   elecs=setdiff(arrayObj.getNeighbors(stimElecs,4),stimElecs);
   %elecs=489
   %Art=NaN*zeros(34,512,55);
   %Art(1:34,ind,:)=Arts;
     for i=30:35
            figure(i)
[Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,Art(:,ind,1:Tmax),x,var0);
Apred=NaN*zeros(1,512,Tmax);
Apred(1,ind,:)=Apred1;

subplot(2,2,3)
plot(squeeze(Art(i,elecs,1:Tmax))'-squeeze(Apred(1,elecs,:))')
title('Kernel Interpolation')
subplot(2,2,1)
plot(squeeze(Art(i,elecs,1:Tmax))'-squeeze(Art(i-1,elecs,1:Tmax))')
title('Naive extrapolation')
axis([0 60 -10 30])

subplot(2,2,2)


[Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,i,Art(1:i-1,ind,1:Tmax),x,var0);

plot(squeeze(Art(i,elecs,1:Tmax))'-squeeze(Apred(1,elecs,:))')
title('Kernel Extrapolation')


axis([0 60 -10 30])

        end
%var0=1000
elecs=setdiff(arrayObj.getNeighbors(stimElecs,3),stimElecs);
elecs=setdiff(elecs,patternNo);
%elecs=489
for i=2:35
    krondiag0=1;
%var0=1000
for k=1:2
krondiag0=kron(krondiag0,dL{k});
end
krondiaginv=(exp(x(end))*krondiag0*Kers{3}(i,i)+var0).^(-1);


figure(35+i)
subplot(2,3,1)
plot([1:55]/20,squeeze(Art(i,elecs,1:Tmax)))
title('Artifact (mean across traces)')

ArtF(i,ind,:)=FilterArtifactLocal(Kers,Art(1:i,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,krondiaginv);
ArtFF=NaN*zeros(1,512,Tmax);
ArtFF(i,ind,:)=ArtF(i,ind,:);
subplot(2,3,2)
plot([1:55]/20,squeeze(ArtFF(i,elecs,1:Tmax)))
title('Filtered Artifact')

subplot(2,3,3)


plot([1:55]/20,squeeze(Art(i,elecs,1:Tmax))'-squeeze(ArtFF(i,elecs,1:Tmax))')
title('Diference')
subplot(2,3,4)

plot(squeeze(TracesAll(i,1,elecs,1:Tmax))')
title('Trace')
ArtF(i,ind,:)=FilterArtifactLocal(Kers,squeeze(TracesAll(1:i,1,ind,1:Tmax)),[x log(var0)],i,ind,Q,Qt,krondiaginv);
ArtFF=NaN*zeros(1,512,Tmax);
ArtFF(i,ind,:)=ArtF(i,ind,:);
subplot(2,3,5)

plot([1:55]/20,squeeze(ArtFF(i,elecs,1:Tmax)))
title('Filtered Trace')
subplot(2,3,6)
plot([1:55]/20,squeeze(TracesAll(i,1,elecs,1:Tmax))'-squeeze(ArtFF(i,elecs,1:Tmax))')
title('Difference')
end
%plot([1:55]/20,squeeze(ArtFF(i,ind,:))'-squeeze(Art(i,ind,1:Tmax))')
%Apred]=ExtrapolateArtifactCondEl(Kers{3},Q{3},Qt{3},dL{3},35,squeeze(Art(1:34,69,8)),x,var0)
%ArtF(i,ind,:)=FilterArtifact(Kers,Art(1:i,ind,1:Tmax),[x log(var0)],i,ind,Q,Qt,dL);


for i=448:length(setdiff([1:512],77))
    for t=1:40
el=ind(i);
[Apred]=ExtrapolateArtifactCondEl(Kers{3},Q{3},Qt{3},dL{3},35,squeeze(Art(1:34,el,t)),x,var0);
Apreds(1,el,t)=Apred;
    end
end


for i=1:50
    x(end)=i
    [Apred1]=ExtrapolateArtifactCond(Kers,Q,Qt,dL,35,Art(1:35-1,ind,1:Tmax),x,var0);
Apredss(i,:,:)=NaN*zeros(1,512,55);
Apredss(i,ind,:)=Apred1;
end