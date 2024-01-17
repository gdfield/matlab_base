function params = MakeStimKernels(params)
%Gonzalo Mena, 3/2016
Difs=params.patternInfo.Difs;
Diags=params.patternInfo.Diags;
listCurrents=params.patternInfo.listCurrents;
Art=params.patternInfo.Art;
patternNo=params.patternInfo.patternNo;
stimElecs=params.patternInfo.stimElecsRel;
var0=params.patternInfo.var0;

DifsStim{1}=Difs{1};
DifsStim{2}=zeros(size(Difs{3},1));
options=params.global.options;
breakpoints=findBreakStimElecs(listCurrents);
onsetC=[];
if(params.bundle.findBundle)
    if(~isnan(params.bundle.onsBundle)&&params.bundle.onsBundle>1)
    onsetC=params.bundle.onsBundle-1;
    
    end
end


breakpoints{1}=sort(unique([0 onsetC breakpoints{1}' size(Art,1)]));
params.patternInfo.breakpoints=breakpoints;

DiagsStim{1}=Diags{1};
DiagsStim{2}=breakpoints{1}';


for k=1:length(breakpoints{1})-1
DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=Difs{3}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1));
if(DiagsStim{2}(k+1)-DiagsStim{2}(k)>1)
DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))/max(max(DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))));
end

end


x01=[3 2 1];
%x01=params.arrayInfo.x(1:3);
x02=-5;
x03=15;

for br=1:length(DiagsStim{2})-1

    inter=[DiagsStim{2}(br)+1:DiagsStim{2}(br+1)];
    nVar=[3 1];
nVarCum=cumsum([0 nVar]);

DifsStimTemp=DifsStim;
DifsStimTemp{2}=DifsStim{2}(inter,inter);
DiagsStimTemp=DiagsStim;
DiagsStimTemp{2}=DiagsStim{2}(br:br+1);
f2 = @(Art,x)logDetKronStimElec(squeeze(Art(inter,stimElecs(1),:)),[x log(var0)],DifsStimTemp,[1:5],[1 2],DiagsStimTemp,nVar);

factor=1;
cont=1;
xSti=-11*ones(1,5);
stiOld=xSti(4);
xStiHist=[];
while(xSti(4)<-10||xSti(4)>5)
g2 = @(x)f2(Art,x);
if(length(inter)>2)

xSti = fminunc(g2,[x01 factor*x02 x03],options);

else
    xSti=[x01 x02 x03];
    
end
evalsf(cont)=g2(xSti);
xStiHist=[xStiHist;xSti];

factor=factor/2;
if(mod(cont,10)==0)
factor=1;
end

cont=cont+1;
if(cont==10)
contind=find(evalsf==min(evalsf));
xSti=xStiHist(contind(1),:);
break
end
end


%x02=xSti(4);
%x03=xSti(5);

for k=1:2
types=[1 2];

[Ker KerD]=evalKernels(DifsStimTemp{k},DiagsStimTemp{k},xSti(nVarCum(k)+1:nVarCum(k+1)),types(k));
[a b]=eig(Ker);
QStim{br}{k}=a';
QtStim{br}{k}=a;
dLStim{br}{k}=diag(b);
KersStim{br}{k}=Ker;
xStim(br,:)=xSti;
end
end

params.patternInfo.KersStim=KersStim;
params.patternInfo.dLStim=dLStim;
params.patternInfo.QtStim=QtStim;
params.patternInfo.QStim=QStim;
params.patternInfo.DifsStim=DifsStim;
params.patternInfo.DiagsStim=DiagsStim;
params.patternInfo.xStim = xStim;