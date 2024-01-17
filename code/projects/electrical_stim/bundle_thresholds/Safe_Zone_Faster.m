%% Step 1: Create rows of adjacency matrix using graph traversal
function [SafeZone_Indices] = Safe_Zone_Faster(dataset, writeDir, DominantBand)
tic;

patternNos = 1:512;
NumElectrodes = length(patternNos);

% Load matrix containing the electrode numbers for the 512-electrode MEA
temp = load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/arrayPositions512.mat'));
positions = temp.positions;
xCoords = positions(:,1);
yCoords = positions(:,2);

load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/bundle_thresholds/Radius1.mat'));
load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/bundle_thresholds/Radius3.mat'));

% Hardcoded: electrode indices for the borders of the array
topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;
borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);

hexCoords = elecHexCoords();
hexArray = elecHexArray(hexCoords);

% Choose the dataset and path to use
%dataset = '/scratch0/karthik/NeuralRecordings/CompleteSets/11920153/data001-data002/';
%DominantBand = 1; % Human knows dominant band (direction against axons)
% AxonDirection = 2; % Only for computing max-hops: 1 is up-down, 2 left-right - will compute on fly

% Create an empty matrix for the adjacency matrix
Adjacency = zeros(length(patternNos), length(patternNos));

% In the process of looking at the data, we will extract the number of
% stimulation amplitudes, number of trials, time samples etc...
NumAmplitudes = 0;
NumSamples = 0;
NumRepeats = 0;

% Go through data from each stimulating electrode.
% Populate the rows of the adjacency matrix using A-star

% Also do a growth test using A-star NOW. This will reduce dimensionality
% of data ALL the upcoming steps, since we will eliminate all stimulation
% amplitudes which fail the growth test.

for i = patternNos
    disp(['Step 1: Extracting Adjacency Matrix Row ' num2str(i)]);
    movieNos = [];
    patternNo = patternNos(i); % current electrode number
    pathToAnalysisData = dataset;
    patternNoString = ['p' num2str(patternNo)];
    files = dir([pathToAnalysisData patternNoString]);
    
    % Create a list of the movie numbers for this electrode
    for j = 1:length(files)
        if strfind(files(j).name, patternNoString) == 1
            mIndices = strfind(files(j).name, 'm');
            movieNos = [movieNos str2double(files(j).name(mIndices(end)+1:end))];
        end
    end
    movieNos = sort(movieNos);
    mIndices = 2:size(movieNos,2);
    
    % Take the data at the lowest amplitude to subtract recording artifact
    % Dimensions: Repeats x Recording Electrodes x Samples
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ...
        patternNo, movieNos(1), 99999);
    [nextA,nextB,nextC] = size(dataTraces);
    
    % Keep these fields updated - they will be needed for CoM calculation
    NumRepeats = max(nextA, NumRepeats);
    NumAmplitudes = max(size(movieNos,2), NumAmplitudes);
    NumSamples = max(nextC, NumSamples);
    
    firstArtifact = mean(dataTraces,1);
    
    % Now go to the last stimulation amplitude
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ...
        patternNo, movieNos(mIndices(end)), 99999);
    [amps, stimChan, ~] = getStimAmps(pathToAnalysisData, ...
        patternNo, movieNos(mIndices(end)));
    
    % Create an artifact matrix that can be subtracted
    subtractionMatrix = repmat(firstArtifact, [size(dataTraces,1) 1]);
    
    % Subtract the recording artifact - also ignore beginning part.
    % Dimensions: Repeats x Recording Electrode x Time
    [meanData, minIndices] = min(mean(dataTraces(:,:,10:40)-subtractionMatrix(:,:,10:40), 1), [], 3);
    
    % This will give you a set of electrodes for the path - use these
    % electrodes and the dataTraces to populate adjacency matrix
    axonPath = findAxonBundlePathFastEdited(meanData, topBorder, rightBorder, bottomBorder, leftBorder, hexCoords, 0, 0);
    Adjacency(patternNo, axonPath) = squeeze(abs(meanData(axonPath)));
end

% We will NEVER use more than 25 samples
NumSamples = min(NumSamples, 25);

%% Step 2: Symmetrize, and create matrix for clustering
% Method of symmetrization
disp('Step 2: Symmetrization');
A = Adjacency*Adjacency' + Adjacency'*Adjacency;

% Hardcoded: electrode indices for the borders of the array
topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;

% Make the sets disjoint
rightBorder = setdiff(rightBorder, topBorder);

bottomBorder = [505:512 8:8:136 129:135];
bottomBorder = setdiff(bottomBorder, union(topBorder, rightBorder));

leftBorder = 129:8:249;
leftBorder = setdiff(leftBorder,union(union(topBorder,rightBorder),bottomBorder));

borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);

%% Step 3: Perform spectral community detection based on modularity, plot outputs
disp('Step 3: Community Detection');

% Now compute degrees and form matrices that can be used for clustering
degrees = sum(A,2);
D = diag(degrees);

% Number of edges as sum
m = sum(degrees)./2;

% Form matrices
L =  D - A; % Laplacian
LN = eye(NumElectrodes) - D^(-0.5)*A*D^(-0.5); % Normalized Laplacian
B = A - (1/(2*m)).*(degrees*degrees'); % Modularity Matrix

partition = ones(NumElectrodes,1);
modularity = (1./(2*m))*sum(sum(B));
deltamod = 1;

% Recursively partition the electrodes based on max eigenvalue
while deltamod > 0
    sets = unique(partition);
    bestdelta = 0;
    bestpartition = partition;
    for i = 1:1:length(sets)
        % Get max eigenvalues of B for this set in partition
        partitionnew = partition;
        [V1,D1] = eigs(B(partitionnew==i,partitionnew==i),1,'lm');
        
        % Update the new partition
        newpart = partitionnew(partitionnew==i);
        newpart(V1<0) = (length(sets)+1).*ones(length(newpart(V1<0)),1);
        partitionnew(partitionnew==i) = newpart;
        
        % Compute the new modularity
        newdelta = (1/(2*m)).*(sum(sum(B(partitionnew==i,partitionnew==i))) +...
            sum(sum(B(partitionnew==length(sets)+1,partitionnew==length(sets)+1))) -...
            sum(sum(B(partition==i,partition==i))));
        
        % If this new partition improves the modularity
        if newdelta > bestdelta
            bestdelta = newdelta;
            bestpartition = partitionnew;
        end
    end
    partition = bestpartition;
    modularity = modularity + bestdelta;
    deltamod = bestdelta;
end

%% Step 4 : Clean out spurious sets in the partition
% (dont cross array, or completely disconnected from rest of partition)
disp('Step 4: Clean out spurious sets');

% Sets with valid entry of 1 are not spurious. Other ones are.
Valid = zeros(length(sets),1);

for k = 1:1:length(sets)
    % Check for atleast two distinct borders that the set hits
    thisPartition = patternNos(partition == sets(k));
    
    checks = ~isempty(intersect(thisPartition,topBorder));
    checks = checks + ~isempty(intersect(thisPartition,rightBorder));
    checks = checks + ~isempty(intersect(thisPartition,bottomBorder));
    checks = checks + ~isempty(intersect(thisPartition,leftBorder));
    
    if checks >= 2 % set hits atleast two distinct borders...
        Valid(k) = 1;
    end
end
GoodSets = sets(Valid == 1);

% For spurious sets, each electrode needs to compute its distance to each of
% the non spurious sets, and reassign itself
for k = 1:1:length(sets)
    % If this set is spurious, go through each electrode and reassign
    if ~Valid(k)
        toFix = patternNos(partition == sets(k));
        for j = 1:1:length(toFix)
            currentElec = toFix(j);
            bestDist = 1000000000000000;
            bestPart = sets(k);
            for q = 1:1:length(GoodSets)
                % Extract all the electrodes from the candidate non-spurious set
                candidateSet = patternNos(partition == GoodSets(q));
                temp = sqrt((positions(currentElec,2)-positions(candidateSet,2)).^2 + ...
                    (positions(currentElec,1)-positions(candidateSet,1)).^2);
                
                % If this is the closest non-spurious set, reassign electrode
                if min(temp) < bestDist
                    bestDist = min(temp);
                    bestPart = GoodSets(q);
                end
            end
            partition(currentElec) = bestPart;
        end
    end
end

%% Step 5 : Connectivity testing and fixes
disp('Step 5: Connectivity constraints');
% Now create a new graph (with edges to neighbors), decompose it into
% connected components, and check that each cluster corresponds to a single
% connected component. Otherwise, reassign to different clusters.
G = zeros(length(patternNos), length(patternNos));
for i = 1:1:length(patternNos)
    newedges = setdiff(Radius1(i,:), 0);
    % Make sure the cluster is the same as electrode i
    newedges = newedges.*[partition(newedges)==partition(i)]';
    
    % take set diff with 0 again so we only add edges to the connected part
    newedges = setdiff(newedges, 0);
    G(i, newedges) = 1;
end

% Now find connected components in the graph
[nc, sz, mm] = networkComponents(G);
% For each of the connected components, there must be a cluster that
% matches EXACTLY. Otherwise need to redo the partition.
sets = unique(partition);

% First fix: If a connected component has length 1, its clearly wrong.
% Reassign ASAP.
for component = 1:1:nc
    currentcheck = mm{component};
    if length(currentcheck) ==1
        setA = Radius1(currentcheck, :);
        partition(currentcheck) = mode(partition(setdiff(setA,0)));
    end
end

% Second fix: any remaining disconnected components have length > 1. Can
% fix now with radius after recomputing G
G = zeros(length(patternNos), length(patternNos));
for i = 1:1:length(patternNos)
    newedges = setdiff(Radius1(i,:), 0);
    % Make sure the cluster is the same as electrode i
    newedges = newedges.*[partition(newedges)==partition(i)]';
    
    % Take set diff with 0 again so we only add edges to the connected part
    newedges = setdiff(newedges, 0);
    G(i, newedges) = 1;
end
[nc, sz, mm] = networkComponents(G);

% For each of the connected components, there must be a cluster that
% matches EXACTLY. Otherwise need to redo the partition
sets = unique(partition);
for component = 1:1:nc % sweep over the connected components
    currentcheck = mm{component}; % Current connected component for checking
    valid = 0; % flag for checking if that component matches perfectly
    % with a partite set
    for partite = 1:1:length(sets) % sweep over the partite sets
        toFix = patternNos(partition==sets(partite)); % electrodes from partite set
        % This means the current connected component is exactly a partition
        if isempty(setdiff(toFix,currentcheck))
            valid = 1;
        end
    end
    
    % Now if valid, we found a matching connected component. Otherwise,
    % create a new set in the partition based on modularity
    if ~valid
        % Check if two borders are hit by the partition
        % Check for at least two distinct borders that the set hits
        thisComponent = patternNos(currentcheck); % Electrodes from component
        checks = ~isempty(intersect(thisComponent,topBorder));
        checks = checks + ~isempty(intersect(thisComponent,rightBorder));
        checks = checks + ~isempty(intersect(thisComponent,bottomBorder));
        checks = checks + ~isempty(intersect(thisComponent,leftBorder));
        
        if checks >= 2
            twoBorders = 1;
        else
            twoBorders = 0;
        end
        
        if twoBorders % If there are two borders hit, we need to just create
            % a new set with the proper index.
            next = setdiff(1:1:max(partition),sets);
            if ~isempty(next)
                partition(currentcheck) = next(1);
            else
                partition(currentcheck) = max(partition)+1;
            end
        else
            % If there are not two borders that are hit, we need to just
            % find the partition which encompasses the two electrodes?
            setA = Radius1(thisComponent,:);
            setB = union(0,thisComponent);
            Outwards = setdiff(setA,setB);
            hitPartitions = unique(partition(Outwards));
            partition(thisComponent) = mode(partition(setdiff(setA,0)));
        end
    else
        % Check for hitting two borders
        thisComponent = patternNos(currentcheck); % electrodes from component
        
        checks = ~isempty(intersect(thisComponent,topBorder));
        checks = checks + ~isempty(intersect(thisComponent,rightBorder));
        checks = checks + ~isempty(intersect(thisComponent,bottomBorder));
        checks = checks + ~isempty(intersect(thisComponent,leftBorder));
        
        if checks >= 2
            twoBorders = 1;
        else
            twoBorders = 0;
        end
        
        if twoBorders
            % If two borders hit, and partition is connected, leave it
        else
            % If there are not two borders that are hit, just
            % find the partition which encompasses the isolated electrodes
            setA = Radius1(thisComponent,:);
            setB = union(0,thisComponent);
            Outwards = setdiff(setA,setB);
            partition(thisComponent) = mode(partition(Outwards));
        end
    end
end

%% Step 6: Extract angles of the bands
disp('Step 6: Extract angles of the bands');

types = unique(partition);
% Clean up the partition indices to be contiguous indices
for i = 1:1:length(partition)
    partition(i) = find(types == partition(i));
end

angles = zeros(max(partition),1); % The "angle" that the partition runs in
perpangles = zeros(max(partition),1); % Angle of the perpendicular line
perpslopes = zeros(max(partition),1); % Slope of the perpendicular line

for i=1:1:max(partition)
    if ~isempty(partition(partition == i))
        P = polyfit(xCoords(partition == i),...
            yCoords(partition == i), 1); % Fit least squares line through electrode
        perpslopes(i) = -1./P(1);
        angles(i) = rad2deg(atan2(P(1),1)); % Always results in acute angle!
        perpangles(i) = rad2deg(atan2(-1./P(1),1));
    end
end

% NEW: compute axon direction, show partitions
anglesTrunc = abs(angles(angles ~=0));
AxonDirs = 2.*[anglesTrunc < 30] + 1.*[anglesTrunc >= 30];
AxonDirection = mode(AxonDirs);

% Need to vary colors before plotting.
colors = {'r';'k';'b';'g';'cyan';'magenta';'yellow';[1,0.6,0];...
    [0,0.75,0.75];[0.75,0.75,0];[0.68,0.47,0];[0.48,0.06,0.89];[0.95,0.87,0.73];...
    [0.2,0.4,0];[0.8,0.8,1];[0.5,0.5,0.5]};
[numColors,t] = size(colors);

types = unique(partition);
figure;
hold on;
for l=1:1:length(types)
    scatter(positions(partition==types(l),1),positions(partition==types(l),2),150,colors{types(l)},'filled');
end
pause(1);

%% Step 7: Compute Bands, Center of Mass, velocities (magnitude and angles)
%  and closest electrodes to each CoM

% Each partite set is divided into two "half"-bands around the stimulating
% electrode.
Bands1 = cell(512);
Bands2 = cell(512);

% Center of Mass in each "Half"-Band - X Coordinate, Y Coordinate stored
% separately.
Centers1X = zeros(NumElectrodes, NumAmplitudes, NumSamples);
Centers1Y = zeros(NumElectrodes, NumAmplitudes, NumSamples);

Centers2X = zeros(NumElectrodes, NumAmplitudes, NumSamples);
Centers2Y = zeros(NumElectrodes, NumAmplitudes, NumSamples);

% Store the index of the closest electrode to each CoM
closestElecInfo1 = cell(NumElectrodes, NumAmplitudes);
closestElecInfo2 = cell(NumElectrodes, NumAmplitudes);

% Velocities in each "Half"-band hold the magnitude and angle of movement
% of the CoM from sample to sample. First entry is magnitude, second is angle.
Velocities1A = zeros(NumElectrodes, NumAmplitudes, NumSamples);
Velocities2A = zeros(NumElectrodes, NumAmplitudes, NumSamples);

% Relative angle of movement - same as angle but converted to be relative
% to the angle of the partite-set
relAngles1 = zeros(size(Velocities1A));
relAngles2 = zeros(size(Velocities2A));
% Displacement of CoM from stimulating electrode
displacement1 = zeros(NumElectrodes,NumAmplitudes,NumSamples);
displacement2 = zeros(NumElectrodes,NumAmplitudes,NumSamples);

% A nonlinear temporal filter applied to waveforms before computing CoM -
% stores the time sample indices of the 5 first order statistics for each waveform thats captured.
NonlinearFilter5 = cell(NumElectrodes, NumAmplitudes);
stimulationLevels = zeros(NumElectrodes, NumAmplitudes);

for stimElec = 1:1:NumElectrodes
    disp(['Step 7: Bands, Centers, Velocities for pattern ' num2str(stimElec) ' of ' num2str(NumElectrodes)]);
    % Relevant electrodes - in same partition, or partition of neighbors
    RelElecs = ismember(partition, unique(partition(setdiff(Radius1(stimElec,:),0))));
    x = xCoords(RelElecs);
    y = yCoords(RelElecs);
    
    xStim = xCoords(stimElec);
    yStim = yCoords(stimElec);
    
    % Compute the y-intercept of the line perpendicular to partition, that
    % runs through the stimulation electrode
    c_int = yStim - perpslopes(partition(stimElec))*xStim;
    m_slope = perpslopes(partition(stimElec));
    relevant = patternNos(RelElecs);
    y_cap = m_slope.*x + c_int;
    
    % Now decide the two "half"-bands. Band1 are the electrodes above the
    % perpendicular line. Band2 are below. Assign electrodes right on the
    % line to the smaller "half"-band.
    if sum([y >= y_cap]) > sum([y < y_cap])
        Band1 = relevant(y > y_cap);
        Band2 = relevant(y <= y_cap);
    else
        Band1 = relevant(y >= y_cap);
        Band2 = relevant(y < y_cap);
    end
    
    % Always remove the stimElec from both "half"-bands. Its corrupted
    % with artifact and could ruin the CoM calculation.
    Band1 = setdiff(Band1, stimElec);
    Band2 = setdiff(Band2, stimElec);
    
    Bands1{stimElec} = Band1;
    Bands2{stimElec} = Band2;
    
    % Combine closest electrode info here too
    startQ = hexCoords(stimElec,1);
    startR = hexCoords(stimElec,2);
    
    % Now need to access the data for this stimElec and use to get CoM
    movieNos = [];
    patternNoString = ['p' num2str(stimElec)];
    files = dir([dataset patternNoString]);
    
    % Create a list of the movie numbers for this electrode
    for j = 1:length(files)
        if strfind(files(j).name, patternNoString) == 1
            mIndices = strfind(files(j).name, 'm');
            movieNos = [movieNos str2double(files(j).name(mIndices(end)+1:end))];
        end
    end
    movieNos = sort(movieNos);
    mIndices = 2:size(movieNos,2);
    
    % Now take the data at the lowest amplitude to subtract recording artifact
    % Dimensions: Repeats x Recording Electrodes x Samples
    dataTraces=NS_ReadPreprocessedData(dataset, '', 0, ...
        stimElec, movieNos(1), 99999);
    [amps, stimChan, ~] = getStimAmps(dataset, ...
        stimElec, movieNos(1));
    
    % Store the Stimulation Amplitude
    stimulationLevels(stimElec, 1) = abs(amps);
    
    % Create a recording artifact matrix that can be subtracted
    firstArtifact = mean(dataTraces,1);
    subtractionMatrix = repmat(firstArtifact, [size(dataTraces,1) 1]);
    
    for stimAmp = 2:1:size(movieNos,2)
        % Extract the data recorded when stimulating with stimAmp
        % and take the mean across repeats
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, ...
            stimElec, movieNos(stimAmp), 99999);
        [amps, stimChan, ~] = getStimAmps(pathToAnalysisData, ...
            stimElec, movieNos(stimAmp));
        
        % Save the current amplitudes in an array
        stimulationLevels(stimElec, stimAmp) = abs(amps);
        
        subtractionMatrix = repmat(firstArtifact, [size(dataTraces,1) 1]); % size changes for stupid reasons
        
        meanData = squeeze(mean(dataTraces-subtractionMatrix, 1)); % Full data after artifact subtraction
        % Dimensions of meanData: Recording Electrode x Time
        
        meanDataClean = meanData(:,1:25); % Dont count stuff after sample 25.
        
        [B,I] = sort(meanDataClean'); % Sort the data in ascending order
        orderedData = B';
        orderingData = I';
        orderingData = orderingData(:,1:5); % Sample indices of first 5 order statistics
        % Create the Nonlinear Filter - Zero out waveforms except for
        % these 5 samples
        NonlinearFilter5{stimElec, stimAmp} = orderingData;
        
        % Now compute CoM in each "Half"-band at each time sample
        for t = 1:(NumSamples+1)
            % If we are past the first time sample, and the past CoM
            % was not the stimulating electrode, include Radius-3
            % electrodes around previous CoM also in the calculation
            % to better track the propagation
            if t > 1 && (Centers1X(stimElec,stimAmp,t-1)~=xStim || Centers1Y(stimElec,stimAmp,t-1)~=yStim)
                % Never include stimElec or 0
                toInclude = setdiff(Radius3(Closest1,:), [0, stimElec]);
                toInclude = union(toInclude, Closest1);
                toInclude = setdiff(toInclude, [0, stimElec]);
                
                % Re-compare to y_cap after including extra points
                x_toInclude = xCoords(toInclude);
                y_toInclude = yCoords(toInclude);
                y_cap = m_slope.*x_toInclude + c_int;
                
                % Assign electrodes right on the line to the smaller "half"-band.
                if sum([y_toInclude >= y_cap]) > sum([y_toInclude < y_cap])
                    Band1_Current = union(toInclude(y_toInclude > y_cap), Band1);
                else
                    Band1_Current = union(toInclude(y_toInclude >= y_cap), Band1);
                end
            else
                % There was no closest electrode to previous CoM, or
                % this is the first sample
                Band1_Current = Band1;
            end
            
            % Fix: Include past CoM Radius 2 also for Band 2...
            if t > 1 && (Centers2X(stimElec,stimAmp,t-1) ~= xStim || Centers2Y(stimElec,stimAmp,t-1) ~= yStim)
                toInclude = setdiff(Radius3(Closest2,:), [0, stimElec]);
                toInclude = union(toInclude, Closest2);
                toInclude = setdiff(toInclude, [0, stimElec]);
                
                % Recalculate y_cap for this radius 2
                x_toInclude = xCoords(toInclude);
                y_toInclude = yCoords(toInclude);
                y_cap = m_slope.*x_toInclude + c_int;
                
                % Assign the equality band to smaller band
                if sum([y_toInclude >= y_cap]) > sum([y_toInclude < y_cap])
                    Band2_Current = union(toInclude(y_toInclude <= y_cap), Band2);
                else
                    Band2_Current = union(toInclude(y_toInclude < y_cap), Band2);
                end
            else
                Band2_Current = Band2;
            end
            
            if ~isempty(Band1_Current)
                meanData_5 = meanData(:,t).*sum(ismember(NonlinearFilter5{stimElec, stimAmp},t), 2);
                data11 = meanData_5(Band1_Current); % Take only current band elecs
                data1 = -1*data11.*[data11<0];
                
                xDistance1 = xCoords(Band1_Current) - xStim;
                yDistance1 = yCoords(Band1_Current) - yStim;
                
                Centers1X(stimElec,stimAmp,t) = sum(data1.*xDistance1)./sum(data1) + xStim;
                Centers1Y(stimElec,stimAmp,t) = sum(data1.*yDistance1)./sum(data1) + yStim;
                if (isnan(Centers1X(stimElec,stimAmp,t)) || isnan(Centers1Y(stimElec,stimAmp,t)))
                    Centers1X(stimElec,stimAmp,t) = xStim;
                    Centers1Y(stimElec,stimAmp,t) = yStim;
                end
                
                diffs = sqrt((Centers1X(stimElec,stimAmp,t)-xCoords).^2 + (Centers1Y(stimElec,stimAmp,t)-yCoords).^2);
                [M,I] = min(diffs);
                Closest1 = I;
                centerQ = hexCoords(Closest1, 1);
                centerR = hexCoords(Closest1, 2);
                
                closestElecInfo1{stimElec,stimAmp} = [closestElecInfo1{stimElec,stimAmp}, elecHexDist(startQ, startR, centerQ, centerR)];
            else
                Centers1X(stimElec,stimAmp,t) = xStim;
                Centers1Y(stimElec,stimAmp,t) = yStim;
                closestElecInfo1{stimElec,stimAmp} = [closestElecInfo1{stimElec,stimAmp}, 0];
                Closest1 = stimElec;
            end
            
            if ~isempty(Band2_Current)
                meanData_5 = meanData(:,t).*sum(ismember(NonlinearFilter5{stimElec, stimAmp},t), 2);
                data22 = meanData_5(Band2_Current);
                data2 = -1*data22.*[data22<0];
                
                xDistance2 = xCoords(Band2_Current) - xStim;
                yDistance2 = yCoords(Band2_Current) - yStim;
                
                Centers2X(stimElec,stimAmp,t) = sum(data2.*xDistance2)./sum(data2) + xStim;
                Centers2Y(stimElec,stimAmp,t) = sum(data2.*yDistance2)./sum(data2) + yStim;
                if (isnan(Centers2X(stimElec,stimAmp,t)) || isnan(Centers2Y(stimElec,stimAmp,t)))
                    Centers2X(stimElec,stimAmp,t) = xStim;
                    Centers2Y(stimElec,stimAmp,t) = yStim;
                end
                
                diffs = sqrt((Centers2X(stimElec,stimAmp,t)-xCoords).^2 + (Centers2Y(stimElec,stimAmp,t)-yCoords).^2);
                [M,I] = min(diffs);
                Closest2 = I;
                centerQ = hexCoords(Closest2, 1);
                centerR = hexCoords(Closest2, 2);
                
                closestElecInfo2{stimElec,stimAmp} = [closestElecInfo2{stimElec,stimAmp},elecHexDist(startQ, startR, centerQ, centerR)];
            else
                Centers2X(stimElec,stimAmp,t) = xStim;
                Centers2Y(stimElec,stimAmp,t) = yStim;
                closestElecInfo2{stimElec,stimAmp} = [closestElecInfo2{stimElec,stimAmp}, 0];
                Closest2 = stimElec;
            end
            
            % Displacement from the stimulating electrode
            displacement1(stimElec,stimAmp,t) = norm([(Centers1X(stimElec,stimAmp,t)-positions(stimElec,1)) ...
                (Centers1Y(stimElec,stimAmp,t)-positions(stimElec,2))]);
            displacement2(stimElec,stimAmp,t) = norm([(Centers2X(stimElec,stimAmp,t)-positions(stimElec,1)) ...
                (Centers2Y(stimElec,stimAmp,t)-positions(stimElec,2))]);
            
            % Now compute velocities for the previous time point
            if t > 1
                if ~isempty(Band1_Current)
                    ang = rad2deg(atan2(squeeze(Centers1Y(stimElec,stimAmp,t)-Centers1Y(stimElec,stimAmp,t-1)),...
                        squeeze(Centers1X(stimElec,stimAmp,t)-Centers1X(stimElec,stimAmp,t-1))));
                    
                    % Store angles in the range [0,360]
                    Velocities1A(stimElec,stimAmp,t-1) = ang.*[ang >= 0] + [ang < 0].*[ang + 360];
                    relAngles1(stimElec,stimAmp,t-1) = mod(Velocities1A(stimElec,stimAmp,t-1) - perpangles(partition(stimElec)),360);
                    
                end
                if ~isempty(Band2_Current)
                    ang = rad2deg(atan2(squeeze(Centers2Y(stimElec,stimAmp,t)-Centers2Y(stimElec,stimAmp,t-1)),...
                        Centers2X(stimElec,stimAmp,t)-Centers2X(stimElec,stimAmp,t-1)));
                    
                    % convert angles to be in range [0,360]
                    Velocities2A(stimElec,stimAmp,t-1) = ang.*[ang >= 0] + [ang < 0].*[ang + 360];
                    relAngles2(stimElec,stimAmp,t-1) = mod(Velocities2A(stimElec,stimAmp,t-1) - perpangles(partition(stimElec)),360);
                end
            end
        end
    end
end

%% Step 8: Time Window Computation
timeWindow1MA = cell(NumElectrodes, NumAmplitudes);
timeWindow2MA = cell(NumElectrodes, NumAmplitudes);

baseWindow1 = cell(NumElectrodes,1);
baseWindow2 = cell(NumElectrodes,1);

for i = 1:NumElectrodes
    % Keep track of moving average of start and end points
    movAvStart1 = 0;
    movAvEnd1 = 0;
    movAvStart2 = 0;
    movAvEnd2 = 0;
    
    % Array of all Starts and Ends
    movAvStartArray1 = [];
    movAvEndArray1 = [];
    movAvStartArray2 = [];
    movAvEndArray2 = [];
    
    % Counter to keep track of how many good timewindows have been found
    stableTW1ctr = 0;
    stableTW2ctr = 0;
    
    % Find the last amplitude to use
    endAmp = NumAmplitudes-1;
    for j = NumAmplitudes-1:-1:2
        if ~isempty(closestElecInfo1{i,j})
            endAmp = j;
            break;
        end
    end
    
    for j = endAmp:-1:2
        % Band 1
        mInc1 = monotoneByDiff(squeeze(closestElecInfo1{i,j})); % check monotonicity
        mInc1 = mInc1(1:25); % Chop after 25th sample
        aInt1 = angleInterval(squeeze(relAngles1(i,j,:)));
        aInt1(aInt1==-1) = 0; % If this angle is between 180
        
        checkASignMInc1 = (aInt1).*(mInc1)'; % combine both
        checkASignMInc1 = checkASignMInc1(1:25); % Now runs of 0, 10, 1
        
        % Now first ignore the 10's in the sequence (these are the
        % ends/starts of runs)
        temp = [checkASignMInc1==1];
        difftemp = diff(temp');
        [~, check1Starts] = find(difftemp~=0);
        check1Starts = check1Starts+1; % To account for diff having 1 value less than original array
        
        check1Starts = [1 check1Starts]; %
        check1Runs = diff([check1Starts (length(checkASignMInc1)+1)]);
        
        goodStarts1 =  [temp(check1Starts) == 1]; % Make sure the checkASignMInc1 is 1 at the start of the run
        check1Starts = check1Starts(goodStarts1); % Now have only the proper starts
        check1Runs = check1Runs(goodStarts1); % Now have only the proper runs
        
        % Now do the check to add in the previous '10' as part of the start
        % and run, if it exists
        for q = 1:1:length(check1Starts)
            if check1Starts(q) > 1
                previousASign = checkASignMInc1(check1Starts(q)-1);
                if (previousASign == 10)
                    check1Starts(q) = check1Starts(q) - 1;
                    check1Runs(q) = check1Runs(q) + 1;
                end
            end
        end
        
        % Possibly merge runs now
        check1Ends = check1Starts + check1Runs - 1;
        toDelete = [];
        extendedStarts1 = [check1Starts 0];
        extendedEnds1 = [0 check1Ends];
        
        diffMerge = [check1Starts 0] - [0 check1Ends] - 1;
        for a = 2:1:(length(diffMerge)-1)
            if diffMerge(a) <= 2 % There is a small gap between two angle runs
                GapToCheck1 = mInc1(extendedEnds1(a):extendedStarts1(a));
                if (sum(GapToCheck1-ones(size(GapToCheck1))) == 0)
                    toDelete = [toDelete a];
                end
            end
        end
        
        % Remove the small gaps between angle runs
        extendedStarts1(toDelete) = [];
        extendedEnds1(toDelete) = [];
        
        % Finish merging
        check1Starts = extendedStarts1(1:end-1);
        check1Ends = extendedEnds1(2:end);
        check1Runs = check1Ends - check1Starts + 1;
        
        startHops1 = squeeze(closestElecInfo1{i,j}(check1Starts));
        endHops1 = squeeze(closestElecInfo1{i,j}(check1Starts+check1Runs-1));
        jumps1 = endHops1 - startHops1;
        
        filter = ([jumps1>0]'&[startHops1<=3]'); % make sure atleast one jump, start hops makes sense
        check1Starts = check1Starts(filter);
        check1Runs = check1Runs(filter);
        if (movAvEnd1==0) && (movAvStart1==0)
            [M,I] = max(check1Runs);
            if M~=0
                timeWindow1MA{i,j} = check1Starts(I):1:(check1Starts(I)+check1Runs(I)-1);
                
                movAvStart1 = timeWindow1MA{i,j}(1);
                movAvEnd1 = timeWindow1MA{i,j}(end);
                movAvStartArray1 = [movAvStart1 movAvStartArray1];
                movAvEndArray1 = [movAvEnd1 movAvEndArray1];
                
                stableTW1ctr = stableTW1ctr + 1;
            end
        else
            % Now just find the time window that has the largest intersection
            % with the stabilized time window. Use that one.
            intTW = zeros(1,length(check1Starts));
            for q = 1:1:length(check1Starts)
                timestoCheck1 = check1Starts(q):1:(check1Starts(q)+check1Runs(q)-1);
                relevanthops1 = squeeze(closestElecInfo1{i,j}(timestoCheck1));
                intTW(q) = length(intersect(round(movAvStart1):1:round(movAvEnd1), timestoCheck1));
            end
            [M,I] = max(intTW);
            if M~=0 % If there is something with a non-trivial intersection, use it for now
                timeWindow1MA{i,j} = check1Starts(I):1:(check1Starts(I)+check1Runs(I)-1);
                
                if (stableTW1ctr < 3) % use only the last 3 "good" amplitudes for estimating TW
                    movAvStartArray1 = [timeWindow1MA{i,j}(1) movAvStartArray1];
                    movAvEndArray1 = [timeWindow1MA{i,j}(end) movAvEndArray1];
                    movAvStart1 = sum(movAvStartArray1)/length(movAvStartArray1);
                    movAvEnd1 = sum(movAvEndArray1)/length(movAvEndArray1);
                    stableTW1ctr = stableTW1ctr + 1;
                end
            end
        end
        
        % "Half"-Band 2
        mInc2 = monotoneByDiff(squeeze(closestElecInfo2{i,j})); % check monotonicity
        mInc2 = mInc2(1:25); % Chop after 25th sample
        aInt2 = angleInterval(squeeze(relAngles2(i,j,:)));
        aInt2(aInt2 == 1) = 0;
        
        checkASignMInc2 = (aInt2).*(mInc2)'; % combine both
        checkASignMInc2 = checkASignMInc2(1:25); % Runs of 0, -10, -1
        
        % Now first ignore the -10's in the sequence (these are the
        % ends/starts of runs)
        temp = [checkASignMInc2==-1];
        difftemp = diff(temp');
        [~, check2Starts] = find(difftemp~=0);
        check2Starts = check2Starts+1; % to account for diff having 1 value less than original array
        
        check2Starts = [1 check2Starts]; %
        check2Runs = diff([check2Starts (length(checkASignMInc2)+1)]);
        
        goodStarts2 = [temp(check2Starts) == 1]; % Make sure the checkASignMInc1 is -1 at the start of the run
        check2Starts = check2Starts(goodStarts2); % Now have only the proper starts
        check2Runs = check2Runs(goodStarts2); % Now have only the proper runs
        
        % Now do the check to add in the previous '-10' as part of the start
        % and run, if it exists
        for q = 1:1:length(check2Starts)
            if check2Starts(q) > 1
                previousASign = checkASignMInc2(check2Starts(q)-1);
                if (previousASign == -10)
                    check2Starts(q) = check2Starts(q) - 1;
                    check2Runs(q) = check2Runs(q) + 1;
                end
            end
        end
        
        % Possibly merge runs now
        check2Ends = check2Starts + check2Runs - 1;
        toDelete = [];
        
        extendedStarts2 = [check2Starts 0];
        extendedEnds2 = [0 check2Ends];
        
        diffMerge = [check2Starts 0] - [0 check2Ends] - 1;
        for a = 2:1:(length(diffMerge)-1)
            if diffMerge(a) <= 2 % There is a small gap between two angle runs
                GapToCheck2 = mInc2(extendedEnds2(a):extendedStarts2(a));
                if (sum(GapToCheck2-ones(size(GapToCheck2))) == 0)
                    toDelete = [toDelete a];
                end
            end
        end
        
        % Remove the small gaps between angle runs
        extendedStarts2(toDelete) = [];
        extendedEnds2(toDelete) = [];
        
        % Finish merging
        check2Starts = extendedStarts2(1:end-1);
        check2Ends = extendedEnds2(2:end);
        check2Runs = check2Ends - check2Starts + 1;
        
        startHops2 = squeeze(closestElecInfo2{i,j}(check2Starts));
        endHops2 = squeeze(closestElecInfo2{i,j}(check2Starts+check2Runs-1));
        jumps2 = endHops2 - startHops2;
        
        filter = ([jumps2>0]'&[startHops2<=3]'); % Need atleast 1 jump, start near stim elec.
        check2Starts = check2Starts(filter);
        check2Runs = check2Runs(filter);
        if (movAvEnd2==0) && (movAvStart2==0)
            [M,I] = max(check2Runs);
            if M~=0
                timeWindow2MA{i,j} = check2Starts(I):1:(check2Starts(I)+check2Runs(I)-1);
                movAvStart2 = timeWindow2MA{i,j}(1);
                movAvEnd2 = timeWindow2MA{i,j}(end);
                movAvStartArray2 = [movAvStart2 movAvStartArray2];
                movAvEndArray2 = [movAvEnd2 movAvEndArray2];
                
                stableTW2ctr = stableTW2ctr + 1;
            end
        else
            % Now just find the time window that has the largest intersection
            % with the stabilized time window. Use that one.
            intTW = zeros(1,length(check2Starts));
            for q = 1:1:length(check2Starts)
                timestoCheck2 = check2Starts(q):1:(check2Starts(q)+check2Runs(q)-1);
                relevanthops2 = squeeze(closestElecInfo2{i,j}(timestoCheck2));
                intTW(q) = length(intersect(round(movAvStart2):1:round(movAvEnd2), timestoCheck2));
            end
            [M,I] = max(intTW);
            if M~=0 % If there is something with a non-trivial intersection, use it for now
                timeWindow2MA{i,j} = check2Starts(I):1:(check2Starts(I)+check2Runs(I)-1);
                if (stableTW2ctr < 3) % use only the last 3 "good" amplitudes for estimating TW
                    movAvStartArray2 = [timeWindow2MA{i,j}(1) movAvStartArray2];
                    movAvEndArray2 = [timeWindow2MA{i,j}(end) movAvEndArray2];
                    movAvStart2 = sum(movAvStartArray2)/length(movAvStartArray2);
                    movAvEnd2 = sum(movAvEndArray2)/length(movAvEndArray2);
                    stableTW2ctr = stableTW2ctr + 1;
                end
            end
        end
    end
    
    baseWindow1{i} = round(movAvStart1):1:round(movAvEnd1);
    baseWindow2{i} = round(movAvStart2):1:round(movAvEnd2);
    disp(['Step 8: Time Windows for pattern ' num2str(i)]);
end

%% Step 9: Feature Extraction

startTimes = zeros(NumElectrodes,NumAmplitudes,2);
endTimes = zeros(NumElectrodes,NumAmplitudes,2);
totalJump = zeros(NumElectrodes,NumAmplitudes,2);
startHop = zeros(NumElectrodes,NumAmplitudes,2);

for i = 1:NumElectrodes
    for j = NumAmplitudes-1:-1:2
        if ~isempty(timeWindow1MA{i,j})
            startTimes(i,j,1) = timeWindow1MA{i,j}(1);
            endTimes(i,j,1) = timeWindow1MA{i,j}(end);
            totalJump(i,j,1) = closestElecInfo1{i,j}(timeWindow1MA{i,j}(end)) - closestElecInfo1{i,j}(timeWindow1MA{i,j}(1));
            startHop(i,j,1) = closestElecInfo1{i,j}(timeWindow1MA{i,j}(1));
        end
        if ~isempty(timeWindow2MA{i,j})
            startTimes(i,j,2) = timeWindow2MA{i,j}(1);
            endTimes(i,j,2) = timeWindow2MA{i,j}(end);
            totalJump(i,j,2) = closestElecInfo2{i,j}(timeWindow2MA{i,j}(end)) - closestElecInfo2{i,j}(timeWindow2MA{i,j}(1));
            startHop(i,j,2) = closestElecInfo2{i,j}(timeWindow2MA{i,j}(1));
        end
    end
    disp(['Step 9: Feature Extraction for pattern ' num2str(i)]);
end

%% Step 10: Bundle detection using CoM

limit_TW1MA_NoWeight = 2*ones(NumElectrodes,1);
limit_TW2MA_NoWeight = 2*ones(NumElectrodes,1);
endAmplitudes = zeros(NumElectrodes,1);

parfor i = 1:NumElectrodes
    disp(['Step 10: Bundle detection B1 for pattern ' num2str(i) ' no weights.']);
    endAmp = NumAmplitudes-1;
    for j = NumAmplitudes-1:-1:2
        if ~isempty(timeWindow1MA{i,j})
            endAmp = j;
            endAmplitudes(i) = j;
            break;
        end
    end
    inter = timeWindow1MA{i,endAmp};
    for j = endAmp:-1:2
        inter = intersect(inter, timeWindow1MA{i,j});
        if (length(inter) <= 1)
            limit_TW1MA_NoWeight(i) = j+1;
            break;
        end
    end
    
    disp(['Step 10: Bundle detection B2 for pattern ' num2str(i) ' no weights.']);
    endAmp = NumAmplitudes-1;
    for j = NumAmplitudes-1:-1:2
        if ~isempty(timeWindow2MA{i,j})
            endAmp = j;
            endAmplitudes(i) = j;
            break;
        end
    end
    inter = timeWindow2MA{i,endAmp};
    for j = endAmp:-1:2
        inter = intersect(inter, timeWindow2MA{i,j});
        if (length(inter) <= 1)
            limit_TW2MA_NoWeight(i) = j+1;
            break;
        end
    end
end

%% Step 12: Computation of the max hops for each electrode (used for jumps)
disp('Step 11: Max Hops Computation');
MaxHops = zeros(NumElectrodes,2);
for e = 1:1:NumElectrodes
    if AxonDirection == 1 % compute max hops in the up-down direction
        I = find(unique(yCoords) == yCoords(e));
        MaxHops(e,1) = length(unique(yCoords)) - I;
        MaxHops(e,2) = I;
    else
        % Fix the row, so look at hops in terms of 32 electrodes
        currentRow = unique(xCoords(yCoords == yCoords(e)));
        I = find(currentRow == xCoords(e));
        if angles(partition(e)) < 0
            MaxHops(e,1) = I;
            MaxHops(e,2) = length(currentRow) - 1;
        else
            MaxHops(e,1) = length(currentRow) - 1;
            MaxHops(e,2) = I;
        end
    end
end

%% Step 12: Growth - we've used A-star to build matrices. Orthogonal test.
SpreadData = cell(NumElectrodes,1);
Growth_res = zeros(NumElectrodes,1);
Growth_IDX = zeros(NumElectrodes,1);
parfor i=1:1:NumElectrodes
    disp(['Step 12: Growth test for pattern ' num2str(i)]);
    SpreadTemp = zeros(endAmplitudes(i),1);
    [bundleMeans, bundleElecs, bundleTimes] = getBundleVoltagesAStarStartEndFromEndAmpl(dataset, ...
        i, 0);
    SpreadTemp(2:endAmplitudes(i)) = squeeze(bundleMeans(1:endAmplitudes(i)-1,1));
    SpreadData{i} = SpreadTemp;
    data = SpreadData{i};
    cummean = cumsum(data')./(1:length(data));
    sq_err = zeros(length(data),1);
    for j = 2:length(data);
        sq_err(j) = (data(j) - cummean(j-1))^2;
    end
    idx = find(sq_err >= 60, 1);
    if isempty(idx) % no bundle activation found
        Growth_res(i) = inf;
        Growth_IDX(i) = inf;
    else
        Growth_res(i) = stimulationLevels(i,idx);
        Growth_IDX(i) = idx;
    end
end
Growth_res = round(100.*Growth_res)./100;

%% Step 13: Logic function for estimation, and jumps

% Bundle thresholds and jump information for each "Half"-band
Band1_res = zeros(NumElectrodes,1);
Band1_idx = zeros(NumElectrodes,1);
Jump1_res = cell(NumElectrodes,1);

Band2_res = zeros(NumElectrodes,1);
Band2_idx = zeros(NumElectrodes,1);
Jump2_res = cell(NumElectrodes,1);

parfor i = 1:1:NumElectrodes
    Band1_res(i) = stimulationLevels(i,limit_TW1MA_NoWeight(i));
    Jump1_res{i} = squeeze(totalJump(i,1:endAmplitudes(i),1));
    % Also keep track of jumps sequence, so can use in test
    Band1_idx(i) = limit_TW1MA_NoWeight(i);
    
    Band2_res(i) = stimulationLevels(i,limit_TW2MA_NoWeight(i));
    Jump2_res{i} = squeeze(totalJump(i,1:endAmplitudes(i),2));
    % Also keep track of jumps sequence, so can use in test
    Band2_idx(i) = limit_TW2MA_NoWeight(i);
end
% Round to 2 decimal places so no dumb false negatives
Band1_res = round(100.*Band1_res)./100;
Band2_res = round(100.*Band2_res)./100;

% Now, based on dominant band decide which side to give more precedence to
if DominantBand == 1
    Dominant_res = Band1_res;
    Nondominant_res = Band2_res;
    
    Dominant_Hops = squeeze(MaxHops(:,1));
    Nondominant_Hops = squeeze(MaxHops(:,2));
    
    Dominant_jumps = Jump1_res;
    Nondominant_jumps = Jump2_res;
    
    Dominant_idx = Band1_idx;
    Nondominant_idx = Band2_idx;
else
    Dominant_res = Band2_res;
    Nondominant_res = Band1_res;
    
    Dominant_Hops = squeeze(MaxHops(:,2));
    Nondominant_Hops = squeeze(MaxHops(:,1));
    
    Dominant_jumps = Jump2_res;
    Nondominant_jumps = Jump1_res;
    
    Dominant_idx = Band2_idx;
    Nondominant_idx = Band1_idx;
end

% Logic Function and Jumps application
Bidir = zeros(NumElectrodes,1); % This will store the amplitudes of bundle thresholds
Bundle_IDX = zeros(NumElectrodes,1); % This will store indices of stim.

parfor i = 1:1:NumElectrodes
    if Dominant_Hops(i) <= 2 % CASE 0 : Too close to border
        Bidir(i) = NaN;
    elseif Dominant_Hops(i) >= Nondominant_Hops(i)
        indexD = JumpsFunction(Dominant_jumps{i}, Dominant_idx(i), Dominant_Hops(i));
        indexND = JumpsFunction(Nondominant_jumps{i}, Nondominant_idx(i), Nondominant_Hops(i));
        if indexND ~= inf && ~isnan(indexND)
            ndlevel = stimulationLevels(i, indexND);
        else
            ndlevel = inf;
        end
        
        if indexD ~= inf && ~isnan(indexD)
            dlevel = stimulationLevels(i, indexD);
        else
            dlevel = inf;
        end
        
        if dlevel <= ndlevel % CASE 1: Use Dominant...
            if dlevel <= Growth_res(i)
                Bidir(i) = dlevel;
                Bundle_IDX(i) = Dominant_idx(i);
            else
                Bidir(i) = NaN; % Use growth test otherwise
                Bundle_IDX(i) = NaN;
            end
        else % CASE 2: Compare indexND, indexD with growth test
            if max(dlevel,ndlevel) <= Growth_res(i)
                Bidir(i) = dlevel;
                Bundle_IDX(i) = Dominant_idx(i);
            else
                Bidir(i) = NaN; % Use growth test otherwise
                Bundle_IDX(i) = NaN;
            end
        end
    else % ND Hops larger
        indexD = JumpsFunction(Dominant_jumps{i}, Dominant_idx(i), Dominant_Hops(i));
        indexND = JumpsFunction(Nondominant_jumps{i}, Nondominant_idx(i), Nondominant_Hops(i));
        if indexND ~= inf && ~isnan(indexND)
            ndlevel = stimulationLevels(i, indexND);
        else
            ndlevel = inf;
        end
        
        if indexD ~= inf && ~isnan(indexD)
            dlevel = stimulationLevels(i, indexD);
        else
            dlevel = inf;
        end
        if dlevel <= ndlevel % CASE 3: use ND
            if ndlevel <= Growth_res(i)
                Bidir(i) = ndlevel;
                Bundle_IDX(i) = Nondominant_idx(i);
            else
                Bidir(i) = NaN;
                Bundle_IDX(i) = NaN;
            end
        else % CASE 2: Compare indexND, indexD with growth test
            if max(dlevel,ndlevel) <= Growth_res(i)
                Bidir(i) = max(dlevel,ndlevel);
                Bundle_IDX(i) = max(Dominant_idx(i),Nondominant_idx(i));
            else
                Bidir(i) = NaN; % Use growth test otherwise
                Bundle_IDX(i) = NaN;
            end
        end
    end
    if Bidir(i) == 0
        Bidir(i) = NaN;
        Bundle_IDX(i) = NaN;
    end
end
% Take minimum with growth test
Bidir = min(Bidir, Growth_res);
% If Inf for index, keep it. Otherwise, the minimum index is correct
Bundle_IDX = min(Bundle_IDX, Growth_IDX);
Bundle_IDX(isinf(Bidir)) = Inf;
Bundle_IDX(isnan(Bidir)) = Inf;

% Round to 2 decimal places so no dumb false negatives
Bidir = round(100.*Bidir)./100;

% Output the values
SafeZone_Amplitudes = Bidir;
SafeZone_Indices = Bundle_IDX;

save([writeDir 'SafeZone_Indices.mat'], 'SafeZone_Indices');
toc
end
