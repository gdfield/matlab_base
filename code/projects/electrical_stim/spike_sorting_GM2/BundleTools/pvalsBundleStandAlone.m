function [pval]=pvalsBundleStandAlone(Art,patternNo,findBundleTimes,nNeighborsExclude,arrayObj)
% Gonzalo Mena, 3/2016
% Lauren Grosberg added arrayObj

[Res]=ResidualsElectrodeSimple(Art,patternNo,findBundleTimes);   
pats=arrayObj.getNeighbors(patternNo,nNeighborsExclude);
nConds=size(Art,1);

pval = nan(nConds,1); 
for k=1:nConds
    [~, p]=kstest((log(Res(k,setdiff([1:512],pats)))-nanmean(log(Res(k,setdiff([1:512],pats)))))./nanstd((log(Res(k,setdiff([1:512],pats))))));   
    pval(k)=log(p);
end
pval(1) = NaN;
