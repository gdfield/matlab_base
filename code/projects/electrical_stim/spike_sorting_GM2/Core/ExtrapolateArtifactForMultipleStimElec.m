function [Apreds]=ExtrapolateArtifactForMultipleStimElec(Art,sigma0,xStim,KersStim,listAmps,Polarity,listAmpsNew,PolarityNew,MaxAmp,stimElec)


c0=0;

if(~isequal(Polarity,PolarityNew))
    Art=-Art;
end
if(isequal(listAmps,listAmpsNew))
    Apred=Art;
    return
end


c0=1;
low=0;
begincond=1;
for l=1:length(KersStim)
    nconds=size(KersStim{l}{2},1);
    listAmpsAux=listAmps(c0:c0+nconds-1);
    high=max(listAmpsAux);

    listAmpsNewAux=listAmpsNew(intersect(find(listAmpsNew>low),find(listAmpsNew<=high)));
    
    
    if(isempty(listAmpsNewAux))
         c0=c0+nconds;
         low=high;
        continue
    end
    
    clear d
    clear d2
    for i=1:length(listAmpsNewAux)
        for j=1:length(listAmpsAux)
            d(i,j)=abs(listAmpsNewAux(i)-listAmpsAux(j));
            
            
        end
    end
    
    
    for i=1:length(listAmpsAux)
        for j=1:length(listAmpsAux)
            d2(i,j)=abs(listAmpsAux(i)-listAmpsAux(j));
            
            
        end
    end
    
    
    d=d/MaxAmp;
    d2=d2/MaxAmp;
    lambda=exp(xStim(l,4));
    KersAux=KersStim{l};
    KersExtra{1}=KersAux{1};
    
    KersExtra{2}=(1+sqrt(3)*(d)*lambda).*exp(-sqrt(3)*(d)*lambda);
    KersAux{2}=  (1+sqrt(3)*(d2)*lambda).*exp(-sqrt(3)*(d2)*lambda);
    
    
    krondiag=1;
    for k=1:2
        
        [a b]=eig(KersAux{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
        
        krondiag=kron(krondiag,dL{k});
        
        
    end
    
    
    krondiag=exp(xStim(l,end))*krondiag+sigma0;
    
    prod1=KronProd(Q,squeeze(Art(c0:c0+nconds-1,stimElec,:)));
    krondiaginv=krondiag.^(-1);
    prod1=krondiaginv.*prod1;
    alpha=KronProd(Qt,prod1);
    
    Apred=exp(xStim(l,end))*KronProd(KersExtra,alpha);
    Apreds(begincond:begincond+length(listAmpsNewAux)-1,:,:)=reshape(Apred,length(listAmpsNewAux),1,size(Art,3));
    
    
    begincond=begincond+length(listAmpsNewAux);
    c0=c0+nconds;
    low=high;
end