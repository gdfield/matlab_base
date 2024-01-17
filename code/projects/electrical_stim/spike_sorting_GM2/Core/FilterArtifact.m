function ArtF=FilterArtifactLocal(Kers,Art,x,cmax,ind,Q,Qt,dL)
%Filter the Artifact (mean of spike substracted traces accross repetitions
%of the same stimulus amplitudes) using the model and noise level. Use for
%all electrodes but the stimulating
%Here, x(end) contains log variance sigma0.
%Gonzalo Mena, 03/2016

var0=exp(x(end));

varsmax=length(x);
logs=0;

Kers{3}=Kers{3}(1:cmax,1:cmax);

krondiag0=1;
%var0=1000
for k=1:3
    if(k==3)
        
        [a, b]=eig(Kers{k});
        Q{k}=a';
        Qt{k}=a;
        dL{k}=diag(b);
    end
    krondiag0=kron(krondiag0,dL{k});
end
krondiaginv=(krondiag0+var0*exp(-x(end-1))).^(-1);


Art=Art(1:cmax,:,:);
%Q{3}=1;
%Qt{3}=1;
Kers2=Kers;
Kers2{3}=Kers{3}(cmax,1:cmax);

prod1=KronProd(Q,Art);
prod1=krondiaginv.*prod1;
alpha=KronProd(Qt,prod1);


Filtered=KronProd(Kers2,alpha);
ArtF=reshape(Filtered,1,length(ind),size(Art,3));
