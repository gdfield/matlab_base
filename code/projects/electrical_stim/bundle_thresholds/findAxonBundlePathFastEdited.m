function path = findAxonBundlePathFastEdited(voltageMap, topEdge, rightEdge, bottomEdge, leftEdge, hexCoords, startElec, endElec)

load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/bundle_thresholds/Radius2.mat'));
load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/bundle_thresholds/Radius3.mat'));
load(fullfile(matlab_code_path,'code/projects/electrical_stim/resources/bundle_thresholds/Neighbors2.mat'));

gscore = zeros(512,1);
fscore = zeros(512,1);
borderElecs = [topEdge, rightEdge, bottomEdge, leftEdge];
hexArray = elecHexArray(hexCoords);

% Find the electrode with the greatest voltage (abs) along the edge
if nargin < 8
    if nargin < 7
        startElec = 0;
    end
    endElec = 0;
end

% Start at the electrode that picked up largest amplitude signal on border
if ~startElec
    borderVoltages = voltageMap(borderElecs);
    [~, ind] = min(borderVoltages); 
    maxElec = borderElecs(ind);
    if ~isempty(intersect(topEdge, maxElec))
        % candidateEndings = setdiff(borderElecs, topEdge);
        borderElecs(ismember(borderElecs, topEdge)) = [];
    end
    if ~isempty(intersect(rightEdge, maxElec))
        % candidateEndings = setdiff(borderElecs, rightEdge);
        borderElecs(ismember(borderElecs, rightEdge)) = [];
    end
    if ~isempty(intersect(bottomEdge, maxElec))
        % candidateEndings = setdiff(borderElecs, bottomEdge);
        borderElecs(ismember(borderElecs, bottomEdge)) = [];
    end
    if ~isempty(intersect(leftEdge, maxElec))
        % candidateEndings = setdiff(borderElecs, leftEdge);
        borderElecs(ismember(borderElecs, leftEdge)) = [];
    end
    startElec = maxElec;
end

% Figure out ending electrode - pick the largest electrode on the border
% without including neighbors up to radius 3?.
if ~endElec
    % Get rid of the neighbors up to certain radius
    excludeNeighbors = union(union(union(startElec, ...
        Neighbors2(startElec,:)), Radius2(startElec,:)), Radius3(startElec,:)); 

    borderElecs(ismember(borderElecs, excludeNeighbors)) = [];
    borderVoltages = voltageMap(borderElecs);
    [~, ind] = min(borderVoltages);
    destElec = borderElecs(ind);
    endElec = destElec;
end

voltageMap = abs(voltageMap);
maxVoltage = max(voltageMap, [], 2);

% weighted map is now positive? 
weightedMap = (-voltageMap + maxVoltage)/maxVoltage;

% Get the ending electrode x and y coordinates
destQ = hexCoords(endElec, 1);
destR = hexCoords(endElec, 2);

visited = zeros(512, 1);
openset = zeros(512, 1);
camefrom = zeros(512, 1);

gscore(startElec) = 0;
fscore(startElec) = gscore(startElec) + elecHexDist(hexCoords(startElec, 1), hexCoords(startElec, 2), destQ, destR)*1.001;
openset(startElec) = 1;

while any(openset)
    current = find(fscore == min(fscore(openset == 1)), 1);
    if current == endElec
        path = reconstructPath(camefrom, endElec);
        break;
    end
    
    openset(current) = 0;
    visited(current) = 1;
    neighbors = hexNeighborsFast(current, hexArray, hexCoords);
    for neighbor = neighbors
        if (visited(neighbor))
            continue;
        end
        gestimate = gscore(current) + weightedMap(neighbor); %elecHexDist(hexCoords(current, 1),hexCoords(current, 2), hexCoords(neighbor, 1), hexCoords(neighbor, 2));
        if ~openset(neighbor) || gestimate < gscore(neighbor)
            camefrom(neighbor) = current;
            gscore(neighbor) = gestimate;
            hueristic = elecHexDist(hexCoords(neighbor, 1),hexCoords(neighbor, 2), destQ, destR);
            fscore(neighbor) = gscore(neighbor) + hueristic*1.001;
            openset(neighbor) = 1;
        end
    end
end