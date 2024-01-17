function [bundleMeans] = getBundleVoltagesAStarStartEndFromEndAmpl_Fake(modifiedData)
% Calculate the bundle voltage for each movie of each pattern in
% using the A* pathfinding algorithm to locate the bundle


% Hardcoded: electrode indices for the borders of the array
topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;

hexCoords = elecHexCoords();
hexArray = elecHexArray(hexCoords);

% load in the fake data
modData = modifiedData.data;
stimElec = modifiedData.stimElec;

% Extract data parameters from the mat file
[~, ~, ~, stimAmpSize] = size(modData);

% Bundlemeans - large enough to hold the info for all amps
bundleMeans = zeros(stimAmpSize, 1);


% Take the data at the lowest amplitude to subtract recording artifact
ignore = [];

% Now take last ampl
dataTraces = squeeze(modData(:, :, :, stimAmpSize));


% Ignore the neighbors of the stimulating electrode?
if isempty(ignore)
    ignore = hexNeighborsFast(stimElec, hexArray, hexCoords);
end


% ignore beginning part.
% Ordering: Repeats x Recording Electrode x Time
[meanData, ~] = min(mean(dataTraces(:,:,10:40),1), [], 3);

axonPath = findAxonBundlePathFastEdited(meanData, topBorder, rightBorder, bottomBorder, leftBorder, hexCoords, 0, 0);
startPath = axonPath(1);
endPath = axonPath(end);

for ii = 2:stimAmpSize % now for each of these files we need to run the algorithm
    dataTraces = squeeze(modData(:, :, :, ii));
   
    % Ignore the neighbors of the stimulating electrode?
    if isempty(ignore)
        ignore = hexNeighborsFast(stimElec, hexArray, hexCoords);
    end
    
    % subtract the recording artifact - also ignore beginning part.
    % Ordering: Repeats x Recording Electrode x Time
    [meanData, minIndices] = min(mean(dataTraces(:,:,10:40), 1), [], 3);
    
    axonPath = findAxonBundlePathFastEdited(meanData, topBorder, rightBorder, bottomBorder, leftBorder, hexCoords, startPath, endPath);
    
    bundle = zeros(size(axonPath, 2), 3);
    
    % Mean deflection along the bundle path
    bundle(:, 1) = meanData(axonPath);
    
    % Electrode indices along the bundle path
    bundle(:, 2) = axonPath;
    
    % Timings for the minimums along the path
    bundle(:, 3) = minIndices(axonPath);
    
    % In the structure, remove the electrodes not on the axon path??
    for j = 1:size(bundle, 1)
        if ismember(bundle(j, 2), ignore)
            bundle(j, 1) = NaN;
        end
    end
    
    % Take the mean deflection along the bundle path?
    bundleMean = mean(bundle(~isnan(bundle(:, 1))), 1);
    
    % Store the bundle mean in the structure?
    bundleMeans(ii-1) = bundleMean;

end
end