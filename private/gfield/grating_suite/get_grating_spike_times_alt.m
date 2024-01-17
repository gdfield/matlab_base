function spike_response = get_grating_spike_times_alt(spike_times, stim_struct, varargin)

p = InputParser;
p.addParameter('accumulate_type', 'direction', @isstring)

p.addParameter('direction', [], @isnumeric);
p.addParameter('TP', [], @isnumeric);
p.addParameter('SP', [], @isnumeric);
p.addParameter('contrast', []);
p.addParameter('background', []);

p.parse(varargin{:});


% extract the number of times each stimulus was presented.
num_repeats = datarun.stimulus.repetitions;
num_stim = length(datarun.stimulus.combinations);


% initialize output cell of spike times
spike_times = cell(num_rgcs, num_repeats);


switch p.Results.accumulate_type
    
    
    case direction
        
        if isempty(p.Results.TP) || isempty(p.Results.SP)
            error('temporal period and spatial period must be specified')
        end
            
        num_dirs = length(stim_struct.params.DIRECTION);
        
        % find the numeric index to the stim matching the user's input
        stim_indices = zeros(num_dirs, 1);
        for g_dir = 1:num_dirs
            for gstim = 1:num_stim
                if stim_struct.combinations(gstim).TEMPORAL_PERIOD == p.Results.TP && ...
                    stim_struct.combinations(gstim).DIRECTION == stim_struct.params.DIRECTION(g_dir) && ...
                    stim_struct.combinations(gstim).SPATIAL_PERIOD == p.Results.SP
                stim_indices(g_dir) = gstim;
                break
             end
        end
        
        