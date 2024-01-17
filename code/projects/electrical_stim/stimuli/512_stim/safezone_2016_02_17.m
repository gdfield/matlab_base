%% Script to generate the safe zone tests for 2016-02-17
  
% Designed to measure the combination of electrode pairs that are
% perpendicular AND parallel to the axon bundle to determine a model for
% safe zone predictions. 
% Test 7 different amplitude ratios. 
% SET TRIGGER INTERVAL TO 0.5 SECONDS

clear all;

% as seen in vision-interactive (EI plot):
% 
% 2   |   1
%-----------
% 3   |   4

%% Define inputs
include_single_elecs = true;
quadrant = 1234;  % 1, 2, 3, 4, 12, or 34 (12 does something in the middle of 1&2, 34 does 3 and 4 together)
pairOrientation = 'downleft'; % horizontal, downleft, downright, vertical (horizontal must be for 60 µm and vertical must be for 30 µm)
secondOrientation = 'downright'; 
delayInMs = 7.5;  % Interval between pulses
saveFiles = 0;    % Set to 1 to save stimulus files, 0 for testing
saveName = 'safezone_quad1234_2016_02_17_9'; %Descriptive name for the stimulus files 
pairRatios = [1    1   1    0.75 0.5 0.25 1;
              0.25 0.5 0.75 1    1   1    1]; % size(pairRatios)=2,number of ratios
%% Define electrodes variable
elec_coords = electrode_positions(512);

if quadrant == 1
    electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)>=0);
elseif quadrant == 2;
    electrodes = find(elec_coords(:,1)<0  & elec_coords(:,2)>=0);
elseif quadrant == 3;
    electrodes = find(elec_coords(:,1)<0 & elec_coords(:,2)<0);
elseif quadrant == 3.1
    electrodes = find(elec_coords(:,1) <40 & elec_coords(:,1)>-940 ...
        & elec_coords(:,2)<0);
elseif quadrant == 34;
    electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
        & elec_coords(:,2)<0);
elseif quadrant == 12;
    electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
        & elec_coords(:,2)>0);
elseif quadrant == 1234;
    electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
        & elec_coords(:,2)>-225 & elec_coords(:,2)<225);
else
    electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)<0);
end

%% Define array variable to produce patterned stimuli.

nSamples=10000; %0.5 s
timeShift=0;
delay=delayInMs*20;

pPerChunk = floor(nSamples/delay);

array = zeros(length(electrodes), 0);
elec_1 = zeros(1,0); %keeps track of electrode 1 (for ordering)
stim_type = zeros(1,0);
for r = 1:size(pairRatios,2)
for ii = 1:length(electrodes)
    rel_coords = zeros(size(elec_coords));
    rel_coords(:,1) = elec_coords(:,1) - elec_coords(electrodes(ii),1);
    rel_coords(:,2) = elec_coords(:,2) - elec_coords(electrodes(ii),2);

    switch lower(pairOrientation)
        case{'horizontal'}
            elec_2 = find(rel_coords(:,1)>55 & rel_coords(:,1)<65 & rel_coords(:,2) == 0);           
        case{'downright'}
            elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>-65 & rel_coords(:,2)<-55);
        case{'downleft'}
            elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>55 & rel_coords(:,2)<65);
        otherwise
            error('invalid electrode pair orientation')
    end
    if isempty(elec_2)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because there is no electrode in chosen direction'])
    elseif ~ismember(elec_2, electrodes)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because other electrode isn''t in this quadrant'])
    else
        array = [array zeros(length(electrodes),1)]; %#ok<AGROW>
        array(ii, end) = pairRatios(1,r);
        array(electrodes==elec_2, end) = pairRatios(2,r);
        elec_1 = [elec_1 electrodes(ii)]; %#ok<AGROW>
        stim_type = [stim_type 1]; %#ok<AGROW>
    end
end
end
% Add second orientation. 
for r = 1:size(pairRatios,2)
for ii = 1:length(electrodes)
    rel_coords = zeros(size(elec_coords));
    rel_coords(:,1) = elec_coords(:,1) - elec_coords(electrodes(ii),1);
    rel_coords(:,2) = elec_coords(:,2) - elec_coords(electrodes(ii),2);

    switch lower(secondOrientation)
        case{'horizontal'}
            elec_2 = find(rel_coords(:,1)>55 & rel_coords(:,1)<65 & rel_coords(:,2) == 0);           
        case{'downright'}
            elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>-65 & rel_coords(:,2)<-55);
        case{'downleft'}
            elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>55 & rel_coords(:,2)<65);
        otherwise
            error('invalid electrode pair orientation')
    end
    if isempty(elec_2)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because there is no electrode in chosen direction'])
    elseif ~ismember(elec_2, electrodes)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because other electrode isn''t in this quadrant'])
    else
        array = [array zeros(length(electrodes),1)]; %#ok<AGROW>
        array(ii, end) = pairRatios(1,r);
        array(electrodes==elec_2, end) = pairRatios(2,r);
        elec_1 = [elec_1 electrodes(ii)]; %#ok<AGROW>
        stim_type = [stim_type 1]; %#ok<AGROW>
    end
end
end

if include_single_elecs
    array = [array eye(length(electrodes))];
    elec_1 = [elec_1 electrodes'];
    stim_type = [stim_type 3*ones(1,length(electrodes))];
end

figure; imagesc(array)
xlabel('Pattern Number'); ylabel('electrode index (not electrode number)');

% temp = load('axonBundleThresholds_2016_02_17_5.mat'); 

temp = load('axonBundleThresholds_2016_02_17_9.mat');
bundleThresholds = temp.axonBundleThresholds_2016_02_17_9; 
meanThresh = mean(bundleThresholds(find(bundleThresholds))); 
 
bundleThresholds(bundleThresholds == 0) = meanThresh; 
if size(bundleThresholds,1) < size(bundleThresholds,2)
    bundleThresholds = bundleThresholds'; 
end
bundleThresholdsM = repmat(bundleThresholds,1,size(array,2));
arrayScaled = array.*bundleThresholdsM;

figure; imagesc(arrayScaled); 
title('Scaled array'); 
xlabel('Pattern Number'); ylabel('electrode index (not electrode number)');

% Randomly order patterns.
pattern_order_all = randperm(size(array,2)); 

times = 0:delay:size(array,2)*(delay-1);
nMovieChunks = ceil(times(end)/nSamples);
allPatterns=[];
patternchunklength = ceil(length(pattern_order_all)/nMovieChunks); 
pPerChunk = floor(nSamples/delay);
Chunks=[];
for i=1:nMovieChunks; 
    Start=(i-1)*patternchunklength;
    try
        PatternsForMovie=pattern_order_all(Start+1:Start+patternchunklength);
        Times = (0:delay:delay*(length(patternchunklength)-1)) + timeShift;
    catch
        PatternsForMovie=pattern_order_all(Start+1:end);
        Times = (0:delay:delay*(length(pattern_order_all(Start+1:end))-1)) + timeShift;
    end
    mChunk=NS_MovieChunkGenerationForExperiment(times,nSamples,PatternsForMovie);
    allPatterns = [allPatterns mChunk]; %#ok<AGROW>
end
MovieChunksFile=[nMovieChunks allPatterns];


%% Save files

if saveFiles
    fid = fopen([saveName '_el'],'wb','ieee-le.l64');
    fwrite(fid,electrodes,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_mv'],'wb','ieee-le.l64');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_pt'],'wb','ieee-le.l64');
    fwrite(fid,arrayScaled,'double');
    fclose(fid);
end

%% Simulate the stimulus

figure; 
for pat = pattern_order_all; 
    cla;
    scatter(elec_coords(:,1),elec_coords(:,2),'k'); axis image; 
    e = find(arrayScaled(:,pat)); 
    hold on; 
    scatter(elec_coords(electrodes(e),1),elec_coords(electrodes(e),2),100,...
        arrayScaled(e,pat),'filled'); axis image;
    caxis([0 2.5]); 
    title(num2str(pat)); 
    pause(0.1); 
end