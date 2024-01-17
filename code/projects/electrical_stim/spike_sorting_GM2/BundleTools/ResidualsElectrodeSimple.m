function [Res]=ResidualsElectrodeSimple(Art,patternNos,times)
%Compute the variances (accross time) for each electrode and condition
%Gonzalo Mena, 03/2016
Art0=Art(:,:,times);

Res=NaN*zeros(size(Art,1),size(Art,2));

for j=1:size(Art,2)
    
    for m=1:size(Art,1)
        
        Res(m,j)=nanvar(squeeze(Art0(m,j,:)));
        
    end
end
Res(:,patternNos)=NaN;