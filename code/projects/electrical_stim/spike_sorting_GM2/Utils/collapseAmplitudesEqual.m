function [listCurrentRangesNew listAmpsNew indexesCollapse] = collapseAmplitudesEqual(listAmpsOld, listCurrentRangesOld)

cont=1;
for m=1:size(listAmpsOld,1)
    if(m==1)
        indexesCollapse{1}=1;
        listAmpsNew(1,:)=listAmpsOld(1,:);
        listCurrentRangesNew(1,:)=listCurrentRangesOld(1,:);
        continue
    end
    
    findA=[];
    for k=1:size(listAmpsOld,2)
        if(isempty(find(listAmpsOld(m,k)==listAmpsNew(:,k))))
            findA= [];
            break
            
        else
            
            findA=[findA find(listAmpsOld(m,k)==listAmpsNew(:,k))'];
        end
    end
    
    findAmp=[];
    findAU=unique(findA);
    for k=1:length(findAU)
        a=find(findAU(k)==findA);
        if(length(a)==size(listAmpsNew,2))
            findAmp=findAU(k);
        end
    end
    
    if(isempty(findAmp))
        cont=cont+1;
        indexesCollapse{cont}=m;
        
        listAmpsNew=[listAmpsNew;listAmpsOld(m,:)];
        listCurrentRangesNew=[listCurrentRangesNew;listCurrentRangesOld(m,:)];
          
          else
          indexesCollapse{cont}=[indexesCollapse{cont} m];
    end
end

