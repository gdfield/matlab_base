% For a datarun with single cone RFs, compute the input to each cone.  
%
% set time_offset to change the time offset between spikes and stimuli.  Usually -1 is best.
%
% S = # stimuli
% C = # cones
% N = # cell ids
%
% Results are stored in three variables:
%
% cone_inputs - S x C matrix, each entry a cone activation value
% spike_times - N x S matrix, each entry number of spikes during that stimulus
% spike_rate  - length N cell, each entry a list of spike times (in units of stimulus number)
%
% cone_inputs are saved to a text file located at "save_path"
%

% test_blur = true
write_flag = false;

% LOAD DATARUN


if ~exist('datarun','var')

    switch 8
        case 8;  datarun = load_data('/Users/gdfield/Analysis/2011-10-25-5/data001/data001');
    end

    % load movie_xml_path & other info
    datarun = load_index(datarun);
    
    % load spikes times and trigger times
    datarun = load_neurons(datarun);
    
    % load java object of movie
    datarun = load_java_movie(datarun); 
    
    % load cone weights matrix
    load([single_cone_path '2011-10-25-5_data001_data001-localmax-V2-st4.5' '/Wc.mat'])
end

% Set path to write saved data
save_path = [single_cone_path,'2011-10-25-5_data001_data001-localmax-V2-st4.5','/cone_stimulus.txt'];



% SET PARAMETERS

start_time = 0;
end_time = 1800;
time_offset = -1;


% refresh time
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;

% compute start and stop times in terms of stimulus number
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
stims = start_stim:end_stim;

% note stimulus size
field_width = datarun.stimulus.java_movie.getWidth;
field_height = datarun.stimulus.java_movie.getHeight;


datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, [single_cone_path '2011-10-25-5_data001_data001-bayes-msf_10.00-BW-2-5']);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = get_sta_summaries(datarun, {1,2,3,4,5}, 'keep_stas', false);


% IDENTIFY SPIKES IN TIME-SHIFTED STIMULUS BINS


% get list of cells
cell_indices = get_cell_indices(datarun,'all');

% initialize storage variables
spike_rate = zeros(length(cell_indices),length(stims));
spike_times = cell(length(cell_indices),1);


% go through each cell
for cc = 1:length(cell_indices)

    cell_index = cell_indices(cc);

    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    spike_rate_ = circshift(spike_rate_(start_stim:end_stim),time_offset);

    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end

    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_;
    spike_times{cc} = sort(spike_times_);
    
end



% COMPUTE INPUTS TO EACH CONE

% show output
fprintf('Computing cone input in frames %d to %d...',start_stim,end_stim)
T=text_waitbar;
start_time_ = clock; % note when it started

% initialize storage variables
cone_inputs = zeros(length(stims),size(Wc,2));

clear gauss_params
gauss_params.center_radius = 5;
gauss_params.center_scale = 1;
gauss_params.normalize = 'sum';
gauss_params.x_size = 50;
gauss_params.y_size = 50;
gauss_params.center = [25 25];
image_filt = make_gaussian(gauss_params);

% cycle through each stimulus
for ss = 1:length(stims)

    T=text_waitbar(T,ss/length(stims) - 0.01);



    % note which stim
    this_stim = stims(ss);

    % get new frame
    STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
   new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;
    %new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]);

%     if test_blur
%         filt_frame = zeros(field_width, field_height, 3);
%         for dm = 1:size(new_frame,3);
%             tmp_frame = conv2(new_frame(:,:,dm), image_filt, 'same');
%             filt_frame(:,:,dm) = tmp_frame;
%         end
%         new_frame = filt_frame;
%     end


    %new_frame = new_frame - 0.5;
    new_frame = reshape(new_frame,[],1);

    % convert to cone space
    cone_inputs(ss,:) = full(Wc'*double(new_frame));

end


% display how long it took
fprintf('done (%0.1f seconds)\n',etime(clock,start_time_));


% save out the data to a text file
if write_flag
    fprintf('writing cone_stimulus.txt...')
    dlmwrite(save_path, cone_inputs, '\t')
    fprintf('complete \n')
end

cd ~/Analysis/
cone_fid = fopen('cone_input.bin','w+', 'b');
fwrite(cone_fid, cone_inputs, 'float32');
fclose(cone_fid);

%write spike times for EJ
% NOTE: you may need to make the 'spike_times' directory
for cc = 1:length(cell_indices)
    temp_cell_id = num2str(datarun.cell_ids(cc));
    tmp_save_path = ['~/Analysis/spike_times/',temp_cell_id,'.txt'];
    dlmwrite(tmp_save_path, spike_times{cc}, '\t')
end
    






