function plotActivation(elecRespAuto,figureNumber)
% PLOTACTIVATION() plots the activation curves of cells analyzed using the
% automatic spike sorting algorithm.
% Inputs:
%    elecRespAuto: output structure of automatic spike sorting that
%                  contains sorting results 
%    figureNumber: specify figure number to draw the plot, 0 opens a new
%                  figure window
% Lauren Grosberg 6/2016

if figureNumber ~=0
    figure(figureNumber); clf; 
else
    figure;
end;
for n = 1:length(elecRespAuto.neuronInfo.spikes)
    nspikes= elecRespAuto.neuronInfo.spikes{n};
    size(nspikes)
    nspikes(nspikes>0) = 1;
    probs = nansum(nspikes,2);
    hold on;  plot(elecRespAuto.stimInfo.listAmps,probs./elecRespAuto.stimInfo.nTrials','.-','MarkerSize',20);
end
legend(num2str(elecRespAuto.neuronInfo.neuronIds))
xlabel('stimulation amplitude');
ylabel('activation probability'); ylim([0 1])
prepIdx = find(elecRespAuto.path.pathToPreparation==filesep,2,'last'); 
prep = elecRespAuto.path.pathToPreparation(prepIdx(1)+1:prepIdx(2)-1); 
drIdx = find(elecRespAuto.path.pathToAnalysisData == filesep,2,'last'); 
datarun = elecRespAuto.path.pathToAnalysisData(drIdx(1)+1:drIdx(2)-1); 
title([prep ' ' datarun ' p' num2str(elecRespAuto.stimInfo.patternNo) ' algorithm output']);
end