function [TracesWhiten templatesWhiten var0 varsE]=WhitenTracesTemplates(TracesAll,templates,params,DifsWhiten,DiagsWhiten,varsE)

for n=1:length(templates)
    Time(n)=size(templates{n},2);
end

for n=1:length(templates)
    templates{n}=templates{n}(:,1:min(Time));
end

if(isempty(varsE))
for e=1:size(TracesAll,3);
   aux=squeeze(TracesAll(1,:,e,:));
   aux=aux-repmat(nanmean(aux),size(aux,1),1);
   for j=1:size(aux,1)
       varsE2(j,:)=nanvar(aux(j,:));
   end
   varsE(e)=nanmedian(varsE2);
end
end

xwhiten=params.arrayInfo.xwhiten;
Tmax=params.global.Tmax;

for k=1:2
[Ker KerD]=evalKernels(DifsWhiten{k},DiagsWhiten{k},[xwhiten(k) -100 -100],1);
[r]=chol(inv(Ker));
Kers{k}=r;
end


if(~isempty(templates))
    Ttemplate=size(templates{1},2);
    
Dif1= zeros(Ttemplate); 
DifsWhitenTemplate=DifsWhiten;
DiagsWhitenTemplate=DiagsWhiten;
for jj=1:Ttemplate
    for ii=1:Ttemplate       
        Dif1(ii,jj)=abs(ii-jj);
    end
end

Dif1=Dif1/Tmax;
DifsWhitenTemplate{1}=Dif1;
DiagsWhitenTemplate{1}=[1:Ttemplate]'/Tmax;
KersTemplate=Kers;
[Ker KerD]=evalKernels(DifsWhitenTemplate{1},DiagsWhitenTemplate{1},[xwhiten(1) -100 -100],1);
[r]=chol(inv(Ker));
KersTemplate{1}=r;


for n=1:length(templates)
    templates{n}=templates{n}./repmat(sqrt(varsE'),1,size(templates{n},2));
    templatesWhiten{n}=reshape(exp(xwhiten(end))^(-1/2)*KronProd(KersTemplate,templates{n}),size(templates{n},1),size(templates{n},2));
   
       
end
else
      templatesWhiten=-1;
end


TracesWhiten=TracesAll;



for i=1:size(TracesAll,1)

    for j=1:size(TracesAll,2)
        TracesAll(i,j,:,:)=TracesAll(i,j,:,:)./reshape(repmat(sqrt(varsE'),1,size(TracesAll,4)),1,1,size(TracesAll,3),size(TracesAll,4));
        TracesWhiten(i,j,:,:)=reshape(exp(xwhiten(end))^(-1/2)*KronProd(Kers,squeeze(TracesAll(i,j,:,:))),size(TracesAll,3),size(TracesAll,4));
        
        
    end
    a=TracesWhiten(i,:,params.patternInfo.indRel,:);
    
    varm(i)=nanvar(a(:));
end
if(length(varm)>=5)
var0=nanmean(varm(1:5));
var0=1;
else
    var0 = -1;
    var0=1;
end
