function [params]=InitializeArrayRobustTest(pathToPreparation,arrayObj,varargin)
% Optional inputs:
%           numElecs: 512 or 519 (default 512)
%        folderNames:
% Gonzalo Mena, 03/2016
% Lauren Grosberg edits to implement array object

% Set default params
folderNames = [];
avoidBundleHyper=0;
whiten=0;
typeVarEstimate=2;
subSampleRate=1;
npatterns=15;
useBrownian=0;
useLocalization=[1 1];
nneighbors=5;
% Read the optional input arguments
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

for jj=1:(nbin/2)
    if ~ischar(varargin{jj*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{jj*2-1})
        case 'foldernames'
            folderNames=varargin{jj*2};
        case 'whiten'
            whiten=varargin{jj*2};
        case 'npatterns'
            npatterns=varargin{jj*2};
        case 'avoidbundlehyper'
            avoidBundleHyper=varargin{jj*2};
        case 'typevarestimate'
            typeVarEstimate=varargin{jj*2};
        case 'subsamplerate'
            subSampleRate=varargin{jj*2};
        case 'usebrownian'
            useBrownian=varargin{jj*2};
        case 'uselocalization'    
            useLocalization=varargin{jj*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Load the electrode positions
positions = arrayObj.getPositions;

% Set default parameters
params.global.Tmax=55;
params.global.subSampleRate=1;
params.global.tarray=[0 [7:32]];
params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on');
params.global.nTrial = 50;
params.global.x0=[2.6766    2.6729    1.5639    2.5233    1.9566  -20.7433    0.3893 32.0274];
%params.global.x0=[9.0723 5.1171 2.9485 3.0434 2.3161 -37.4725 0.3176  127.8172];
params.global.positions=positions;
params.global.maxIter=5;
params.global.thresEI=30;
params.global.sortData=0;
params.global.avoidBundleHyper=avoidBundleHyper;
params.global.whiten=whiten;
params.global.saveArt=0;
x0=params.global.x0;
params.global.useStimElec = 0;
params.bundle.findBundle=1;
params.bundle.cutBundle = 0;
params.bundle.nVec=1;
params.bundle.updateFreq =3;
params.bundle.useBundleAlg = 0;
params.bundle.detectThreshold=-12;

params.bundle.findBundleTimes=setdiff(1:params.global.Tmax,[6 7 8]);
params.bundle.nNeighborsExclude=3;

tarray=params.global.tarray;
Tmax=params.global.Tmax;
options=params.global.options;

pathToAnalysisData=[];
patternNo=[];

if isempty(folderNames)
    dirs=dir(pathToPreparation);
    cont=1;
    
    for ii=1:length(dirs)
        if(length(dirs(ii).name)>=4)
            aux=find(dirs(ii).name(1:4)=='data',1);
            if(~isempty(aux))
                folderNames{cont}=dirs(ii).name;
                cont=cont+1;
            end
        end
    end
end

patternList=[];
for f=1:length(folderNames)
    pathAux=fullfile(pathToPreparation, folderNames{f});
    dirs=dir(pathAux);
    
    for ii=3:length(dirs)
        
        aux=find(dirs(ii).name(1)=='p',1);
        if(~isempty(aux))
            if(~isequal(dirs(ii).name(2),'a'))
                
                patternList=[patternList str2num(dirs(ii).name(2:end))];
                
            end
        end
    end
end






%stimElecs=1 %if want to make it work with p2

%shift stimulating elecrode to match the new ordering inducing after
%erasing inactive channels



% Get active electrode list from the array object
activeElecs = arrayObj.getElectrodes;
numActiveElecs = length(activeElecs);
positions=getPositions(arrayObj);
DifPos=arrayObj.difPos;



pattern0 = arrayObj.center;

[theta,rho]  = cart2pol(positions(:,1)-positions(pattern0,1),positions(:,2)-positions(pattern0,2));


els=activeElecs;
rhoUnique=rho(els);
thetaUnique=theta(els);

for f=1:length(rhoUnique)
    for l=1:50
        ArtAll{f,l}=NaN*zeros(Tmax,50);
        cont(l,f)=1;
    end
end

listAmpsAll=[];
var00=[];

sizes=[];
contTraces=1;
while(max(cont)<npatterns)
    patternNo=patternList(unidrnd(length(patternList)));
    
    for i=1:length(folderNames)
        try
            pathToAnalysisData=[pathToPreparation folderNames{i}];
            [TracesAll, Art, var0, listAmps, listCurrents, stimElecs, onset, onsetC, pval, Res, sampledTrials, Trace0]=loadTracesArtSort(pathToAnalysisData,patternNo,Tmax,params.global.nTrial,subSampleRate,arrayObj,'params',params);
            
            if(avoidBundleHyper)
                endCond=onsetC-1;
            else
                endCond=size(Art,1);
            end
            Traces0(contTraces,:,:)=Trace0;
            
            contTraces=contTraces+1;;
            var00=[var00 var0];
            listAmpsAll=[listAmpsAll;listAmps'];
            sizes=[sizes length(listAmps)];
            
            [thetaAux,rhoAux] = cart2pol(arrayObj.xc-arrayObj.xc(patternNo),arrayObj.yc-arrayObj.yc(patternNo));
            thetaAux=thetaAux(activeElecs);
            rhoAux=rhoAux(activeElecs);
            
            for j=1:length(activeElecs)
                
                f=intersect(find(rhoAux(j)==rhoUnique),find(thetaAux(j)==thetaUnique));
                
                if(~isempty(f))
                    
                    if(isnan(endCond)||endCond<=1)
                        for l=1:size(Art,1)
                            ArtAll{f,l}(:,cont(l,f))=[squeeze(Art(l,j,1:Tmax))];
                            cont(l,f)=cont(l,f)+1;
                        end
                    else
                        for l=1:endCond
                            ArtAll{f,l}(:,cont(l,f))=[squeeze(Art(l,j,1:Tmax))];
                            cont(l,f)=cont(l,f)+1;
                        end
                    end
                end
            end
        end
    end
    
end


elsbad=[];
maxCond=find(max(cont')==1)-1;
maxCond=maxCond(1);
for i=1:maxCond
    elsbad=union(elsbad,find(cont(i,:)==1));
end
elsbad=union(els(elsbad),pattern0);

els=setdiff(els,elsbad);

els=intersect(els,arrayObj.getNeighbors(pattern0,nneighbors));

for l=1:maxCond
    for f=1:length(els)
        
        f2=intersect(find(rho(els(f))==rhoUnique),find(theta(els(f))==thetaUnique));
        

        Arts(l,f,:)=nanmean(ArtAll{f2,l}');

    end
end

Dif1= zeros(Tmax,Tmax);
for jj=1:Tmax
    for ii=1:Tmax
        Dif1(ii,jj)=abs(ii-jj);
    end
end
Dif1=Dif1/Tmax;
Difs{1}=Dif1;


Difs{2}=DifPos(els,els)/arrayObj.maxR;


listAmpsAll=unique(listAmpsAll);
params.arrayInfo.listAmpsAll=listAmpsAll;
params.arrayInfo.rho=rho;
Dif3 = zeros(length(listAmpsAll),length(listAmpsAll));
if(useBrownian==0)
    for jj=1:length(listAmpsAll)
        for ii=1:length(listAmpsAll)
            Dif3(ii,jj)=abs(listAmpsAll(ii)-listAmpsAll(jj));
        end
    end
else
    for jj=1:length(listAmpsAll)
        for ii=1:length(listAmpsAll)
            Dif3(ii,jj)=min(listAmpsAll(ii),listAmpsAll(jj));
        end
    end
end


params.arrayInfo.Dif2=DifPos;
params.arrayInfo.Dif3=Dif3;

Dif3=Dif3(1:maxCond,1:maxCond)/max(listAmpsAll);
Difs{3}=Dif3;

Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho(els)/arrayObj.maxR;
Diags{3}=listAmpsAll(1:maxCond)/max(listAmpsAll);

%% Whiten Terms
DifsNoise{1}=Difs{1};

for e=1:numActiveElecs
    maxs(e)=max(max(abs(squeeze(Traces0(:,activeElecs(e),:)))));
end


elsWhiten=activeElecs(find(maxs<median(maxs)));
DifsNoise{2}=DifPos(elsWhiten,elsWhiten)/arrayObj.maxR;
nTrials=contTraces-1;
DifsNoise{3}=eye(nTrials);

DiagsNoise{1}=ones(Tmax,1);
DiagsNoise{2}=ones(length(elsWhiten),1);
DiagsNoise{3}=ones(nTrials,1);

Traces00=Traces0(:,elsWhiten,:);

for e=length(elsWhiten)
    aux=Traces00(:,e,:);
    Traces00(:,e,:)=Traces00(:,e,:)./mean(sqrt(mean(var(squeeze(aux)'))));
end

if(whiten)
    f1=@(Traces0,x)logDetKron(Traces00(:,:,end-Tmax+1:end),[x(1) -100 -100 x(2) -100 -100 0 x(3) -1000],DifsNoise,[1 4 8],[1 1 7],DiagsNoise,[3 3 1]);
    f1=@(Traces0,x)logDetKron(Traces00(:,:,end-Tmax+1:end),[x(1) -100 -100 x(2) -100 -100 0 0 -1000],DifsNoise,[1 4],[1 1 7],DiagsNoise,[3 3 1]);
    
    g1=@(x)f1(Traces0,x);
    xwhiten = fminunc(g1,0*ones(1,2),options);
    xwhiten = [xwhiten 0];
else
    xwhiten=[Inf Inf 0];
end

params.arrayInfo.xwhiten = xwhiten;
params.arrayInfo.DiagsNoise=DiagsNoise;
params.arrayInfo.DifsNoise=DifsNoise;
params.patternInfo.Difs=Difs;
params.patternInfo.Diags=Diags;
params.arrayInfo.Traces0=Traces0;
params.arrayInfo.elsWhiten=elsWhiten;




for i1=1:size(Arts,1);
    
    a=Arts(i1,:,:);
    vars(i1)=nanvar(a(:));
end

if(typeVarEstimate==1)
    var0=mean(var00)/length(var00);
    var0=0;
else
    var0=nanmean(vars(1:5));
    var0=1;
end



if(useBrownian)
    nvar=[3 3 1];
    type=[1 1 8];
  
    if(useLocalization(1)==1)
        if(useLocalization(2)==1)
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1:6) 0  x(end) log(var0)],Difs,setdiff([1:8],[7]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 0 0 -2 0 0 10],options);
            x=[x1(1:6) 0 x1(1)];
        else
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1:4) -100 -100  0 x(end) log(var0)],Difs,setdiff([1:8],[5 6 7]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 0 0 -2 10],options);
             x=[x1(1:4) -100 -100  0 x1(end)];
        end
        
    else
        if(useLocalization(2)==1)
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1) -100 -100 x(2:4) 0  x(end) log(var0)],Difs,setdiff([1:8],[2 3 7]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 -2 0 0 10],options);
            x=[x1(1) -100 -100 x1(2:4) 0 x1(end)];
        else
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1) -100 -100 x(2) -100 -100  0 x(end) log(var0)],Difs,setdiff([1:8],[2 3 5 6 7]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[0 0 10],options);
            x=[x1(1) -100 -100 x1(2) -100 -100  0 x1(end)];
        end
    end
else
    nvar=[3 3 3];
    type=[1 1 1];
  
    if(useLocalization(1)==1)
        if(useLocalization(2)==1)
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1:7) -100 -100  x(end) log(var0)],Difs,setdiff([1:10],[8 9]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 0 0 -2 0 0 -2 10],options);
            x =[x1(1:7) -100 -100  x1(end)];
        else
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1:4) -100 -100 x(5) -100 -100 x(end) log(var0)],Difs,setdiff([1:10],[5 6 8 9]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 0 0 -2 -2 10],options);
            x=[x1(1:4) -100 -100 x1(5) -100 -100 x1(end)];
        end
        
    else
        if(useLocalization(2)==1)
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1) -100 -100 x(2:5) -100 -100  x(end) log(var0)],Difs,setdiff([1:10],[2 3 8 9]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 -2 0 0 -2 10],options);
            x=[x1(1) -100 -100 x1(2:5) -100 -100  x1(end)];
        else
            f1=@(Arts,x)logDetKron(Arts(:,:,:),[x(1) -100 -100 x(2) -100 -100 x(3) -100 -100 x(end) log(var0)],Difs,setdiff([1:10],[2 3 5 6 8 9]),type,Diags,nvar);
            g1=@(x)f1(Arts,x);
            x1 = fminunc(g1,[-2 -2 -2 10],options);
            x=[x1(1) -100 -100 x1(2) -100 -100 x1(3) -100 -100 x1(end)];
            
        end
    end
end

value=g1(x);
params.arrayInfo.value=value;
    
 
  
params.arrayInfo.x= x;
params.arrayInfo.els=els;
params.patternInfo.var0=var0;
params.arrayInfo.patternNo =patternNo;
if(useBrownian==0)
    factp=[0 3 6 9];
    type=[1 1 1];
else
    factp=[0 3 6 7];
    type=[1 1 8];
end

for k=1:3
    [Ker, KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),type(k));
    [a, b]=eig(Ker);
    Q{k}=a';
    Qt{k}=a;
    dL{k}=diag(b);
    Kers{k}=Ker;
end

params.patternInfo.Kers=Kers;
params.patternInfo.Arts=Arts;
