% Script to generate a single electrode scan on the 30um stimulation board.

% Define stimulation time parameters
timeShiftInMs = 7.5;
interPulseLatencyInMs = 7.5;
timeShift = timeShiftInMs*20;
interPulseLatency = interPulseLatencyInMs*20;
numberOfSamples = 10000;

% Hack - the software still assumes that the 512-stim board is the only one
% that exists for electrical stimulation. Therefore, the electrodes file
% (and the Array matrix) can only contain 512 electrodes. 
electrodes = 1:512;
array = eye(512,512); 

% Pattern numbers correspond to the 519 electrode map. 
patterns = 1:519;

% Since 519>512, some channels are disconnected
disconnected=[1 130 259 260 389 390 519];
patternsToUse = setdiff(patterns,disconnected);

% at this point we have the electrodes we want to use for the scan; this
% electrode sequence includes all the electrodes that are actually 
% connected to the electronics (there are 512 of them).
% Please redefine this sequence up to your taste.

% Randomly permute the electrodes in a single chip. Produces a new order
% every time, may be better to define a specific sequence. 
randomChipOrder = randperm(64); 
patternTest = reshape(patternsToUse,64,8)'; 
patternTest(:,randomChipOrder) = patternTest; 
patternsOrder = patternTest(:); 

% Visualize stimulation order to verify it is as you want it
[xc,yc]=getElectrodeCoords519();
figure; scatter(xc(patternsToUse),yc(patternsToUse),100,patternsOrder,'filled'); axis image; axis off; colormap jet; colorbar; 
title('order of electrode scan - complete'); 

figure;
for mc = 1:8
    scatter(xc(patternsOrder((1:64)+64*(mc-1))),yc(patternsOrder((1:64)+64*(mc-1))),100,(1:64),'filled'); axis image; axis off; colormap jet; colorbar;
    title(sprintf('order of electrode scan - movie chunk %0.0f',mc));
    pause(0.5);
end

% Set the sequence order
patternsToUse = patternsOrder; 

% Apply the electrode transform to map the desired 519 patterns on the 512 system
patternsTransformed=NS512_519Inverse(patternsToUse);

%% Generate the stimulation files

times = timeShift:interPulseLatency:timeShift+63*interPulseLatency; % Hard coded to have 64 patterns per movie chunk, can change this if needed.
chunks=[];
p=[];
for ii=1:8
    patternsForMovie = patternsTransformed((ii-1)*64+1:ii*64);
    chunk = NS_MovieChunkGenerationForExperiment(times,numberOfSamples,patternsForMovie);
    chunks =[chunks chunk]; %#ok<AGROW>
end

movieChunksFile=[ii chunks];

keyboard; % Switch to the desired directory for saving the files. 

fid = fopen('519_512_el','wb');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('519_512_pt','wb','ieee-le.l64');
fwrite(fid,array,'double');
fclose(fid);

fid = fopen('519_512_mv','wb','ieee-le.l64');
fwrite(fid,movieChunksFile,'int32');
fclose(fid); 