function datapaths = get_nerve_crush_data(varargin)

% Defines the datasets in the Chichilnisky Lab Nerve Crush Project


p = inputParser;
p.addParameter('cohort', 'control', @(x) ismember(x, {'control', '4day', '7day', '11day'}));
p.parse(varargin{:})

switch p.Results.cohort

    case 'control'

        pathlist{1} = '/Volumes/FieldLab/Chichilnisky/Analysis/controls/20130401A/chunk1/kilosort25/kilosort25';

    case '4day'

        pathlist{1} = [];
        warning('need to fill in these data paths')

    case '7day'

        pathlist{1} = [];
        warning('need to fill in these data paths')

    case '11day'
        
        pathlist{1} = [];
        warning('need to fill in these data paths')

end

datapaths.paths = pathlist;
