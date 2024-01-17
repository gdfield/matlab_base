function [bundleMeans, bundleElecs, bundleTimes] = getBundleVoltagesAStarStartEndFromEndAmpl(path, patternNos, display)
% Calculate the bundle voltage for each movie of each pattern in 
% patternNos, using the A* pathfinding algorithm to locate the bundle

bundleElecs = cell(size(patternNos, 1),50);
bundleTimes = cell(size(patternNos, 1),50);
bundleDeflects = cell(size(patternNos, 1),50);

% Load matrix containing the electrode numbers for the 512-electrode MEA
temp = load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/arrayPositions512.mat'));
positions = temp.positions;

% Hardcoded: electrode indices for the borders of the array
topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;
borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);

hexCoords = elecHexCoords();
hexArray = elecHexArray(hexCoords);

% Bundlemeans why 50? - just large enough to hold the info for all ampls
bundleMeans = zeros(50, 3, size(patternNos, 1));

% Find movie indices for this pattern
for patternIndex = 1:size(patternNos, 2)
    movieNos = [];
    patternNo = patternNos(patternIndex);
    pathToAnalysisData = path; 
    patternNoString = ['p' num2str(patternNos(patternIndex))];
    files = dir([pathToAnalysisData patternNoString]);
 
    for i = 1:length(files)
        if strfind(files(i).name, patternNoString) == 1
            mIndices = strfind(files(i).name, 'm');
            movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))];
        end
    end
    movieNos = sort(movieNos);
    mIndices = 2:size(movieNos,2);
    
    % Order : repeats, recording elctrodes, samples
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ... 
        patternNo, movieNos(1), 99999);
    
    % Take the data at the lowest amplitude to subtract recording artifact
    firstArtifact = mean(dataTraces,1);
    if display
        f = figure; set(f,'Position',[100 360 1000 550]);
        set(f,'Color','white');
    end
    ignore = [];
    
    % Now sweep over all stimulation amplitudes
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ...
        patternNo, movieNos(mIndices(end)), 99999);
    % get stimulus amplitude
    [amps, stimChan, ~] = getStimAmps(pathToAnalysisData, ...
        patternNo, movieNos(mIndices(end)));
    
    % Create an artifact matrix that can be subtracted
    subtractionMatrix = repmat(firstArtifact, [size(dataTraces,1) 1]);
    
    % Ignore the neighbors of the stimulating electrode?
    if isempty(ignore)
        ignore = hexNeighborsFast(stimChan, hexArray, hexCoords);
    end
    % disp(ignore);
    
    % subtract the recording artifact - also ignore beginning part.
    % Ordering: Repeats x Recording Electrode x Time
    [meanData, minIndices] = min(mean(dataTraces(:,:,10:40)-subtractionMatrix(:,:,10:40),1), [], 3);
    
    axonPath = findAxonBundlePathFastEdited(meanData, topBorder, rightBorder, bottomBorder, leftBorder, hexCoords, 0, 0);
    startPath = axonPath(1);
    endPath = axonPath(end);
   
    for  movieIndex = mIndices
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ... 
            patternNo, movieNos(movieIndex), 99999);
        % get stimulus amplitude
        [amps, stimChan, ~] = getStimAmps(pathToAnalysisData, ...
            patternNo, movieNos(movieIndex));
        
        % Create an artifact matrix that can be subtracted
        subtractionMatrix = repmat(firstArtifact, [size(dataTraces,1) 1]);
        
        % Ignore the neighbors of the stimulating electrode?
        if isempty(ignore)
            ignore = hexNeighborsFast(stimChan, hexArray, hexCoords);
        end
        %disp(ignore);
        
        % subtract the recording artifact - also ignore beginning part.
        % Ordering: Repeats x Recording Electrode x Time
        [meanData, minIndices] = min(mean(dataTraces(:,:,10:40)-subtractionMatrix(:,:,10:40), 1), [], 3);
        
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
        bundleMeans(movieIndex-1, 1, patternIndex) = bundleMean;
        
        % Store the stimulation amplitude for this movie
        bundleMeans(movieIndex-1, 2, patternIndex) = amps(1);
        
        % Store the movieNo in the structure as well
        bundleMeans(movieIndex-1, 3, patternIndex) = movieNos(movieIndex);
        
        % Store the Axon Path in another structure
        bundleElecs{patternIndex,movieIndex-1} = axonPath;
        bundleDeflects{patternIndex,movieIndex-1} = axonPath;
        % Store timings of the minimum values along axon in another struct
        bundleTimes{patternIndex,movieIndex-1} = bundle(:, 3);
        
        if display
            cla;
            scatter(positions(:,1), positions(:,2), 350, meanData, 'filled');
            axis off; axis image; colorbar;
            caxis([-50 -10]);
            
            title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)), 'Color', 'black');
            hold on; scatter(positions(bundle(:, 2),1),positions(bundle(:, 2),2),350,[0.5 0.5 0.5], 'filled');
            % text(positions(stimChan,1),positions(stimChan,2),'stimulating electrode')
            pause(0.001);
        end
        
    end
        
end
end