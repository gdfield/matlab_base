function perf = get_pred_fr_metrics( psth, preds, varargin )

%psth : 1 x nbins of actual fr
%preds: number of different predictions x nbins of predicted fr
%optional: 1 or 0 for each metric type [corrcoef rms ev r^2]. default is to compute all of them.

%perf: resulting [corrcoef rms ev r^2]

%future: make it work for many psths at a time, add option for explainable variance

if length(varargin)==0
    mets = [1 1 1 1];
else
    mets = varargin{1};
end

for pr = 1:size(preds,1)
    if mets(1) == 1
        temp = corrcoef(psth,preds(pr,:));
        perf(pr,1) = temp(1,2);
    else
        perf(pr,1) = nan;
    end
    
    if mets(2) == 1
        temp1 = (psth - preds(pr,:)).^2; %diff squared
        perf(pr,2) = sqrt( mean(temp1)); %sqrt of avg
    else
        perf(pr,2) = nan;
    end
    
    if mets(3) == 1
       temp2 = (psth - preds(pr,:)).^2; %diff squared
       temp3 = mean(temp2) / var(psth); %avg / psth var
       perf(pr,3) = 1-temp3; %exp var
    else
        perf(pr,3) = nan;
    end
    
    if mets(4) == 1
        temp4 = mean((psth-preds(pr,:)).^2); %diff squared
        temp5 = mean((psth-mean(psth)).^2); %squared error of psth
        perf(pr,4) = 1-temp4/temp5;
    else
        perf(pr,4) = nan;
    end
        
end




end