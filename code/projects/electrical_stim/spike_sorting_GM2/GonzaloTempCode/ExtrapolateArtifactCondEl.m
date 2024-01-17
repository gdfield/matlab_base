function [Apred]=ExtrapolateArtifactCondEl(Kers,Q,Qt,dL,cmax,Art,x,sigma0)
%Extrapolate artifact at condition cmax given (full) observations until
%condition cmax-1. Extrapolation is done in al electrodes but the
%stimulating
%Gonzalo Mena, 03/2016


%extrapolate
ctest=cmax;
test=ctest;

train=setdiff([1:size(Art,1)],ctest);




Ktrain=Kers(train,train);
Ktraintest=Kers(test,train);


krondiag=1;

        [a b]=eig(Ktrain);
        Q=a';
        Qt=a;
        dL=diag(b);
    
    
    krondiag=kron(krondiag,dL);
    
    



krondiag=krondiag+sigma0*exp(-x(end));
prod1=Q*Art(train);
%prod1=KronProd(Q,Art(train,:,:));
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;
alpha=Qt*prod1;

plot((Ktrain^-1)*Art(train)-alpha)
Apred=Ktraintest*alpha;

Apred=reshape(Apred,length(ctest),size(Art,2),size(Art,3));

