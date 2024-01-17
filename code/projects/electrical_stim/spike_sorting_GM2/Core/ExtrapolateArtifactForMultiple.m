function [Apred]=ExtrapolateArtifactForMultiple(Art,sigma0,x,Kers,indRel,listAmps,Polarity,listAmpsNew,PolarityNew,MaxAmp)
%Extrapolate artifact at condition cmax given (full) observations until
%condition cmax-1. Extrapolation is done in al electrodes but the
%stimulating
%Gonzalo Mena, 03/2016



if(~isequal(Polarity,PolarityNew))
    Art=-Art;
end
if(isequal(listAmps,listAmpsNew))
    Apred=Art;
    return
end





for i=1:length(listAmpsNew)
    for j=1:length(listAmps)
     d(i,j)=abs(listAmpsNew(i)-listAmps(j));
     
    
    end
end


for i=1:length(listAmps)
    for j=1:length(listAmps)
     d2(i,j)=abs(listAmps(i)-listAmps(j));
     
    
    end
end


d=d/MaxAmp;
d2=d2/MaxAmp;
lambda=exp(x(7));
KersExtra{1}=Kers{1};
KersExtra{2}=Kers{2};
KersExtra{3}=(1+sqrt(3)*(d)*lambda).*exp(-sqrt(3)*(d)*lambda);
Kers{3}=  (1+sqrt(3)*(d2)*lambda).*exp(-sqrt(3)*(d2)*lambda); 


 krondiag=1;
for k=1:3

[a b]=eig(Kers{k});
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);

krondiag=kron(krondiag,dL{k});


end


krondiag=exp(x(end))*krondiag+sigma0;

prod1=KronProd(Q,Art(:,indRel,:));
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);

Apred=exp(x(end))*KronProd(KersExtra,alpha);

Apred=reshape(Apred,length(listAmpsNew),length(indRel),size(Art,3));



