function [KersWhitenSRInv]=WhitenKernels(DifsWhiten,DiagsWhiten,xwhiten,nCond)

for k=1:2
[Ker KerD]=evalKernels(DifsWhiten{k},DiagsWhiten{k},[xwhiten(k) -100 -100],1);
[r]=chol(inv(Ker));
KersWhitenSRInv{k}=r;
end

KersWhitenSRInv{3}=eye(nCond);