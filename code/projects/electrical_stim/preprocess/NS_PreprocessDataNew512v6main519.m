function [PDChunkNumber,movieBegin,repetitionNumber,repetitionPeriod,clusterFileName] = NS_PreprocessDataNew512v6main519(pathToRawData, writePath, writePathFigs, arrayID, figureProperties, NS_GlobalConstants, makePlots, patternsToRead, channelsToRead, movieNumberOffset)
% NS_PreprocessDataNew512v6main519(...) preprocesses the raw electrical 
% stimulation data using the 519/30um neurostim board, SB512-3. This
% preprocessing function accounts for an electrode transform that was
% provided to the LabVIEW 512 stimulation software to map the 519/30um
% electrode array map onto the 512/60um electrode array map. This electrode
% remapping occurs 1) before saving the Status files 2) before saving the
% Pattern files and 3) for separating the data into folders labeled by
% stimulation pattern. 
%
% NS_PreprocessDataNew512v6main519(...) also generates and plots the 
% overlaid responses for all patterns and all movies FOR 1-EL. STIMULATION
% (LGrosberg note - did not test this functionality as of 6/2016)
% It is assumed that the number of pattern called in the movie is the 
% number of stimulating electrode.
%
% % % Input arguments:
% pathToRawData, full path to the raw data (string), e.g. '/Volumes/Data/2016-06-13-0/data001';
% writePath, full path to location to save the preprocessed data (string),
%               e.g. '/Volumes/Analysis/2016-06-13-0/data001/';
% writePathFigs, full path to location to save the output figs (string)
% arrayID, array ID used to collect the data. arrayID = 1502 for SB512-3
% figureProperties, e.g.: struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
% NS_GlobalConstants, array-specific structure, contains SamplingFrequency, ChipAddresses, NumberOfChannelsPerChip, CurrentRanges
% makePlots, 0/1, choose whether or not to save plots
% patternsToRead, if movie includes 64 patterns, but we want to process
%                   only some of them, the indexes should be specifiec here. The aray of
%                   indexes is not array of patterns really, but indexes for the array ofPatternsUsed
%                   patterns that will be loaded by this function from the stimulus files. 
% channelsToRead, 
% movieNumberOffset


% % % Output arguments: 
% PDChunkNumber,
% movieBegin,
% repetitionNumber,
% repetitionPeriod,
% clusterFileName

% this function was originally written by Pawel, but has been modified to
% apply the electrode map transformation to the Status and Pattern files by
% Lauren Grosberg, 6/2016. 


% Read in the raw file and the electrode map for the array used
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(pathToRawData);
electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(arrayID);

% Define the umber of samples to save after each stimulus pattern
traceLength = 100;

% Find the paths to LabVIEW output files that contain stimulation info
[parentstr,datarun,~] = fileparts(pathToRawData); 
datarunNumber = datarun(5:end); 
altFileOrg = fullfile(parentstr,'Electrical', 'Output');
if ~isempty(dir(fullfile(parentstr,'movie*')))
    movieFileName= fullfile(parentstr,['movie' datarunNumber]);
    SPfilename    = fullfile(parentstr,['pattern' datarunNumber]);
elseif ~isempty(dir(fullfile(altFileOrg,'movie*')))
    movieFileName= fullfile(altFileOrg, ['movie' datarunNumber]);
    SPfilename    = fullfile(altFileOrg,['pattern' datarunNumber]);
else
    % User selects directory with the labview output files
    correctFileOrg = uigetdir(parentstr, 'Select directory containing LabVIEW output files ''movie*'' and ''pattern*'''); 
    if correctFileOrg
        movieFileName = fullfile(correctFileOrg,['movie' datarunNumber]);
        SPfilename     = fullfile(correctFileOrg,['pattern' datarunNumber]);
    else
        err = MException('MATLAB:missingLabVIEWFiles', ...
            'Cannot find labview output files that define the times and patterns of electrical stimulation');
        throw(err);
    end
end

% Open the movie file, read the header, return the number of movies (number
% of movie chunks x number of amplitudes)
fid0=fopen(movieFileName,'r','b');
header=readMHchunk(fid0);
nMovies = header.number_of_movies

channels=channelsToRead;

nPatternsMax = 0;
nRepetitionsMax = 0;

if arrayID == 1502
    % Load the mapping of the 519/30um stimulation electrodes onto the 512/60um mapping
    correctedElectrodes = get512_519_inverseElecs();
end

% Initialize a stopwatch
t0 = clock;

for ii= 1:1:nMovies-1 % minus 1 because the last movie is expected to be empty. This loop repeats once for each amplitude (in case of scan data).    
    %a ) estimate how much time is left to complete the loop  
    if ii>1 
        finished = (ii-1)/(nMovies-1); % proportion of files created so far
        fprintf('finished writing %0.1f%% of files\n', finished*100);
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        fprintf('estimated time left: %0.1f minutes\n',estimatedTimeLeft/60)
    end

    %b) read in single movie_data_chunk    
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC (command) chunk...
        error('command chunk found in the movie file');
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD (movie data) chunk
        ChunkSize=fread(fid0,1,'int64');
        ChunkData=fread(fid0,ChunkSize,'int32');
        %reading in the movie parameters:
        ChunkData=NS_MovieData(datarunNumber,ii,NS_GlobalConstants,'fullPath',movieFileName);
        [PDChunkNumber,movieBegin,RepetNumber0,repetitionPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);        
        repetitionNumber=min(RepetNumber0,25);
    end    

    %calculate number of events related to specified patterns
    PatternsUsed=MovieData(2:3:length(MovieData));
    
    PatternsIndexesToRead=[]; % 
    for i1=patternsToRead
        PatternIndexes=find(PatternsUsed==i1);  %number of events corresponding to given pattern
        PatternsIndexesToRead=[PatternsIndexesToRead' PatternIndexes']'; 
    end
              
    NumberOfEvents=length(PatternsIndexesToRead);  % this will not work when some patterns are used more than once in the movie!!!
    Events=zeros(NumberOfEvents,1);
    
    %c) read in corresponding pattern_data_chunk, save status into file, predefine array for RAW data  
  
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<ASGLU>
    disp(['ReadPatternDataChunk done for Chunk' num2str(PDChunkNumber) ' ' SPfilename]);
    FullName_status=fullfile(writePath,'status_files', ['status_m' num2str(ii+movieNumberOffset)]);
    if ~exist(fileparts(FullName_status),'dir')
        mkdir(fileparts(FullName_status))
    end
    if arrayID == 1502
        
        NumberOfChips=length(NS_GlobalConstants.ChipAddresses);
        ChipStatus=struct('StimActive',0','ConnectDelay',100,'RefreshDelay',100,'RecordDelay',110,'DischargeDelay',105,'HoldDelay',115);
        ChipsStatus(1:NumberOfChips)=ChipStatus;
        ChannelStatus=struct('active',0,'mode',0,'range',2);
        ChannelsStatus(1:electrodeMap.getNumberOfElectrodes-1)=ChannelStatus;
        remappedStatus=struct('ChipsStatus',ChipsStatus,'ChannelsStatus',ChannelsStatus);
        remappedStatus.ChannelsStatus(correctedElectrodes) = Status.ChannelsStatus; 
        
        % All the values in the status variable are defaults, except for
        % the value in ChannelsStatus.range, which is later used to
        % determine the stimulus amplitude for a particular electrode.
        Status = remappedStatus; 
    end
    save(FullName_status,'Status');
            
    NumberOfPatterns=length(PatternsIndexes);
    nPatternsMax=max(nPatternsMax,NumberOfPatterns);
    length(channelsToRead);
    numberofbytes = NumberOfEvents*repetitionNumber*length(channelsToRead)*traceLength*2 ;   
    Traces=[];
    if numberofbytes>0
        if numberofbytes>100000000 % poniewa? maksymalna wielkosc tablicy to 943 MB, wi?c dla duzych tablic trzeab stworzyc najpierw ma?a tablic?, skonwertowa? na typ int16 i dopiero z takiej sk?ada? du?? talbic? int16
            TracesSmall=int16(zeros(NumberOfEvents,1,length(channelsToRead),traceLength)); 
            Traces=TracesSmall;
            for ble=2:repetitionNumber               
                Traces=[Traces TracesSmall];                
            end
        else
            Traces=int16(zeros(NumberOfEvents,repetitionNumber,length(channelsToRead),traceLength));  
        end                                
    size(Traces);
    %d) read in RAW data for each repetition, for each stimulation event
    %identify pattern number, and save RAW data coresponding to each event into
    %Traces array. We assume each repetition of the movie includes identical
    %collection of patterns (logical, isn't it).

    for j=1:repetitionNumber % one iteration of this loop takes less than 0.5s for 64-channel data (122 events)                
        %RepIndex=j
        RawData=int16(rawFile.getData(movieBegin+repetitionPeriod*(j-1),repetitionPeriod+traceLength)');                        
        for k=1:NumberOfEvents % later on, single iteration should read and organize data for all ... of given pattern in this movie at once (PH, 2010-08-26)
            PatternNumber=PatternsUsed(PatternsIndexesToRead(k));
            t=MovieData(1+(PatternsIndexesToRead(k)-1)*3)  ;                                          
            Events(k)=PatternNumber;
            size(RawData);
            size(RawData(channels+1,t+1:t+traceLength));            
            Traces(k,j,:,:)=RawData(channels+1,t+1:t+traceLength); %this takes less than 0.1ms for 64-channel data                                 
        end        
    end    
    
    DifferentPatterns=unique(Events);
        
    clear RawData;
    
    for l=1:length(DifferentPatterns) %1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        ThisPattern=DifferentPatterns(l);  
        WhichEvents=find(Events==ThisPattern); %which events in this movie corresponded to given pattern
        nRepetitionsMax=max(nRepetitionsMax,repetitionNumber*length(WhichEvents));
        if ~isempty(WhichEvents)            
            Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,ThisPattern);                                    
                                    
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(channels),:),repetitionNumber*length(WhichEvents),length(channels),traceLength);              
            STTS=size(TracesToSave);  
            if arrayID==1502

                TracesToSave519=zeros(STTS(1),519,STTS(3));
                TracesToSave519(:,correctedElectrodes,:)=TracesToSave;
                clear TracesToSave;
                TracesToSave=TracesToSave519;
                clear TracesToSave519;
                STTS=size(TracesToSave);
            end            
        
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);      
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(channels))=channels';
            b(4+length(channels))=length(WhichEvents);
            o=[b' a'];            
                        
            if arrayID==1502               
                correctedPatternNo = correctedElectrodes(ThisPattern);
                FullName=fullfile(writePath,['p' num2str(correctedPatternNo)],['p' num2str(correctedPatternNo) '_m' num2str(ii+movieNumberOffset)]);
            else
                FullName=[writePath filesep 'p' num2str(ThisPattern) filesep 'p' num2str(ThisPattern) '_m' num2str(ii+movieNumberOffset)];
            end
            if ~exist(fileparts(FullName),'dir')
                mkdir(fileparts(FullName));
            end
            fid=fopen(FullName,'wb','ieee-le');                                 
            fwrite(fid,o,'int16');
            fclose(fid);                    
                        
            if arrayID==1502
                FullName_pattern= fullfile(writePath, 'pattern_files',['pattern' num2str(correctedPatternNo) '_m' num2str(ii+movieNumberOffset)]);
                Pattern.channel = correctedElectrodes(Pattern.channel); 
            else
                FullName_pattern=[writePath filesep 'pattern' num2str(ThisPattern) '_m' num2str(ii+movieNumberOffset)];  
            end   
            if ~exist(fileparts(FullName_pattern),'dir')
                mkdir(fileparts(FullName_pattern))
            end
            save(FullName_pattern,'Pattern');
            
            % part below - only for 1-el. scan!!
            if makePlots==1
                ChannelsToPlot=electrodeMap.getAdjacentsTo(ThisPattern,1);                
                TracesToShow=TracesToSave(:,ChannelsToPlot,:);
                %c=int16(mean(TracesToShow));
                %for ci=1:STTS(1)
                %    TracesToShow(ci,:,:)=TracesToShow(ci,:,:)-c;
                %end
            
                y = NS_PlotManySignaturesOnArrayLayoutNew(TracesToShow,ChannelsToPlot,ones(1,STTS(1)),arrayID,figureProperties,NS_GlobalConstants);
                hj=gcf;              
                set(hj,'PaperUnits','inches');
                set(hj,'PaperSize',[13 9]);
                set(hj,'PaperPosition',[0 0 13 9]);
                FullName=[writePathFigs '\' 'p' num2str(ThisPattern) '_m' num2str(ii)];
                print(hj, '-dtiff', '-r120', FullName);
            end
        end
    end                                                                          
    clear Traces;
    clear TracesToSave;
    end
end  
fclose(fid0);
% Create cluster file:
clusterFileName=NS_CreateClusterFile(writePath,datarunNumber,nMovies,nPatternsMax,nRepetitionsMax);