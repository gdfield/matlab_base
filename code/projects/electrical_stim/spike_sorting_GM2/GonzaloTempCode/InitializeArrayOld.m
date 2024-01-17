function params=InitializeArray(pathToPreparation,subSampleRate,arrayObj,varargin)
% Optional inputs:
%           numElecs: 512 or 519 (default 512)
%        folderNames: 
% Gonzalo Mena, 03/2016
% Lauren Grosberg edits to implement array object 

% Set default params
folderNames = []; 
patternNo0 =[]; 

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
        case 'patternno0'
            patternNo0 = varargin{jj*2};
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
params.global.tarray=[0 [7:32]];
params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on');
params.global.nTrial = 50;
params.global.x0=[2.6766    2.6729    1.5639    2.5233    1.9566  -20.7433    0.3893 32.0274];
params.global.positions=positions; 
params.global.maxIter=5;
params.global.thresEI=30;
params.global.sortData=0;
x0=params.global.x0;
params.global.useStimElec = 0;

params.bundle.cutBundle = 0;
params.bundle.nVec=1;
params.bundle.updateFreq =3;
params.bundle.useBundleAlg = 0;
params.bundle.detectThreshold=-12;
params.bundle.findBundle = 0;
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

for f=1:length(folderNames)
    pathAux=fullfile(pathToPreparation, folderNames{f});
    dirs=dir(pathAux);
      
    for ii=3:length(dirs)
        
        aux=find(dirs(ii).name(1)=='p',1);
        if(~isempty(aux))
            if(~isequal(dirs(ii).name(2),'a'))
                if isempty(patternNo0)
                    patternNo=str2num(dirs(ii).name(2:end));
                    pathToAnalysisData=pathAux;
                    break
                else
                    patternNo=str2num(dirs(ii).name(2:end));
                    if(patternNo0==patternNo);
                        pathToAnalysisData=pathAux;
                        break
                    end
                end
            end
        end
    end
end

if(isempty(pathToAnalysisData))
    disp('Could not find useful pattern')
    params.Log{1}='Could not find useful pattern';
    return
end


[TracesAll, Art, var0, listAmps, listCurrents, stimElecs]=loadTracesArtSort(pathToAnalysisData,patternNo0,Tmax,params.global.nTrial,subSampleRate,arrayObj,'params',params);


%stimElecs=1 %if want to make it work with p2

%shift stimulating elecrode to match the new ordering inducing after
%erasing inactive channels
stimElecs=stimElecs-length(find(arrayObj.nullChannels<stimElecs));

if(length(stimElecs)>1)
    disp('only one stimulating electrode supported')
    return
end

Dif1= zeros(Tmax,Tmax); 
for jj=1:Tmax
    for ii=1:Tmax       
        Dif1(ii,jj)=abs(ii-jj);
    end
end
Dif1=Dif1/(max(max(Dif1)));
Difs{1}=Dif1;

% Get active electrode list from the array object
activeElecs = arrayObj.getElectrodes; 
numActiveElecs = length(activeElecs); 

if size(Art,2) ~= numActiveElecs
    disp('size of the artifact should match the number of active channels');
end

Dif2 = zeros(numActiveElecs,numActiveElecs); 
for ii=1:numActiveElecs
    for jj=1:numActiveElecs        
        Dif2(ii,jj)=norm(positions(activeElecs(ii),:)-positions(activeElecs(jj),:),2);
    end
end
positions=positions(activeElecs,:);

[~,rho] = cart2pol(positions(:,1)-positions(stimElecs,1),positions(:,2)-positions(stimElecs,2));

ind=setdiff(1:length(activeElecs),find(rho==0));
rho=rho(ind);

Dif2=Dif2(ind,ind)/max(max(Dif2(ind,ind)));
Difs{2}=Dif2;;

params.global.Dif1 = Dif1;
params.global.Dif2 = Dif2;


Dif3 = zeros(length(listAmps),length(listAmps)); 
for jj=1:length(listAmps)
    for ii=1:length(listAmps)
        Dif3(ii,jj)=abs(listAmps(ii)-listAmps(jj));
    end
end

Dif3=Dif3/max(max(Dif3));
Difs{3}=Dif3;

Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/max(rho);
Diags{3}=listAmps'/max(listAmps);


f1=@(Art,x)logDetKron(Art(:,ind,:),[x(1:7) -100 -100 x(8) log(var0)],Difs,setdiff([1:10],[8 9]),[1 1 1],Diags,[3 3 3]);
g1=@(x)f1(Art,x);

x1 = fminunc(g1,x0,options);
x=[x1(1:7) -100 -100 x1(end)];

params.arrayInfo.x= x;
params.arrayInfo.patternNo =patternNo;
params.arrayInfo.stimElecs= stimElecs;
params.patternInfo.Difs=Difs;
params.patternInfo.Diags=Diags;
params.patternInfo.var0=var0;