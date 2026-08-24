function saverasters(datarun, output_dir, cell_type, varargin)
%SAVERASTERS  Save a consolidated per-cell raster matrix (all trials,
%all matched stimulus conditions) as one CSV file per cell.
%
%   saverasters(datarun, output_dir, cell_type)
%   saverasters(datarun, output_dir, cell_type, 'cell_selection', [31 46 91])
%   saverasters(datarun, output_dir, cell_type, 'direction', 90)
%   saverasters(datarun, output_dir, cell_type, 'direction', 90, 'temporal_period', 256)
%   saverasters(datarun, output_dir, cell_type, 'bin_width', 0.001, 'stim_dur', 8)
%
% REQUIRED INPUTS
%   datarun     - the datarun struct (after load_data/load_neurons/
%                 load_params/load_stim)
%   output_dir  - folder to save the per-cell CSVs into (created if it
%                 doesn't exist)
%   cell_type   - which cell(s) to save. Accepts:
%                   - a string name, matched case-insensitively against
%                     datarun.cell_types{i}.name, e.g. 'ON parasol'
%                   - 'all' to include every cell regardless of type
%                   - a cell array of type-name strings, e.g.
%                     {'ON parasol', 'OFF parasol'} (union of both types)
%                   - a numeric vector of actual cell IDs, e.g. a
%                     mapped_cell_ids list from map.ei - this bypasses
%                     cell-type lookup entirely and uses those cell IDs
%                     directly
%
% OPTIONAL NAME-VALUE INPUTS
%   'cell_selection'    - 'all' (default) to save every cell of the given
%                         cell_type, or a numeric vector of specific cell
%                         IDs (which will be restricted to those also
%                         belonging to cell_type)
%   'direction'         - only include stimulus combinations with this
%                         DIRECTION. Omit to include all directions (or
%                         if your stimulus set has no DIRECTION field).
%   'temporal_period'   - only include combinations with this
%                         TEMPORAL_PERIOD. Omit to include all temporal
%                         periods (or if your stimulus set has no
%                         TEMPORAL_PERIOD field).
%   'bin_width'         - time bin width in seconds for the spike-count
%                         matrix (default 0.001, i.e. 1 ms bins)
%   'stim_dur'          - stimulus duration in seconds (default 8)
%
% BEHAVIOR
%   direction and temporal_period are OPTIONAL and are only applied as
%   filters if your datarun.stimulus.combinations actually has those
%   fields. This lets the same function work on stimulus classes that
%   don't carry direction/temporal-period information at all (e.g. full
%   -field flashes, light steps) - in that case conditions are simply
%   labeled by their combination index instead.
%
%   Each cell's CSV has one row per time bin and one column per
%   (condition, trial) pair, concatenated side by side. Column headers
%   identify the condition and trial, e.g. dir90_TP256_trial1 when
%   direction/TP are available, or cond3_trial1 when they are not.

%% parse inputs
p = inputParser;
addParameter(p, 'cell_selection', 'all');
addParameter(p, 'direction', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'temporal_period', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'bin_width', 0.001, @isnumeric);
addParameter(p, 'stim_dur', 8, @isnumeric);
parse(p, varargin{:});

cell_selection  = p.Results.cell_selection;
direction       = p.Results.direction;
temporal_period = p.Results.temporal_period;
bin_width       = p.Results.bin_width;
stim_dur        = p.Results.stim_dur;

if ischar(cell_selection) || isstring(cell_selection)
    if strcmpi(cell_selection, 'all')
        requested_ids = datarun.cell_ids;
    else
        error('saverasters:badSelection', ...
            'cell_selection must be ''all'' or a numeric vector of cell IDs.');
    end
else
    requested_ids = cell_selection;
end

% resolve cell_type into the set of cell IDs belonging to that type
% accepts: a string name ('ON parasol' or 'all'), a cell array of string
% names ({'ON parasol','OFF parasol'}), or a numeric vector (treated as
% indices into datarun.cell_types, e.g. [1 3])
if ischar(cell_type) || isstring(cell_type)
    if strcmpi(cell_type, 'all')
        type_ids = datarun.cell_ids;
    else
        if ~isfield(datarun, 'cell_types')
            error('saverasters:noCellTypes', ...
                'datarun has no cell_types field - cannot filter by cell_type. Pass ''all'' to skip type filtering.');
        end

        type_names = cellfun(@(x) x.name, datarun.cell_types, 'UniformOutput', false);
        type_match = find(strcmpi(type_names, cell_type));

        if isempty(type_match)
            error('saverasters:badCellType', ...
                'No cell type named "%s" found in datarun.cell_types. Available types: %s', ...
                cell_type, strjoin(type_names, ', '));
        end

        type_ids = datarun.cell_types{type_match(1)}.cell_ids;
    end

elseif iscell(cell_type)
    % cell array of type-name strings -> union of their cell IDs
    if ~isfield(datarun, 'cell_types')
        error('saverasters:noCellTypes', ...
            'datarun has no cell_types field - cannot filter by cell_type. Pass ''all'' to skip type filtering.');
    end

    type_names = cellfun(@(x) x.name, datarun.cell_types, 'UniformOutput', false);
    type_ids = [];

    for tt = 1:length(cell_type)
        this_match = find(strcmpi(type_names, cell_type{tt}));
        if isempty(this_match)
            error('saverasters:badCellType', ...
                'No cell type named "%s" found in datarun.cell_types. Available types: %s', ...
                cell_type{tt}, strjoin(type_names, ', '));
        end
        type_ids = union(type_ids, datarun.cell_types{this_match(1)}.cell_ids);
    end

elseif isnumeric(cell_type)
    % numeric vector -> treated as actual cell IDs directly (e.g. a
    % mapped_cell_ids list from map.ei), not type indices. No cell-type
    % lookup happens in this case; the vector itself IS the cell list.
    type_ids = cell_type(:)';

else
    error('saverasters:badCellTypeArg', ...
        ['cell_type must be a string (e.g. ''ON parasol'' or ''all''), a cell array of ' ...
         'strings (e.g. {''ON parasol'',''OFF parasol''}), or a numeric vector of cell ' ...
         'type indices (e.g. [1 3]).']);
end

% final cell list: requested IDs restricted to the given cell_type
cell_id_list = intersect(requested_ids, type_ids, 'stable');

if isempty(cell_id_list)
    error('saverasters:noCellsMatched', ...
        'No cells matched both cell_type "%s" and the given cell_selection.', cell_type);
end

cell_selection_is_all = (ischar(cell_selection) || isstring(cell_selection)) && strcmpi(cell_selection, 'all');
if ~cell_selection_is_all && length(cell_id_list) < length(requested_ids)
    warning('saverasters:someCellsDropped', ...
        '%d of the requested cell IDs were not of type "%s" and were excluded.', ...
        length(requested_ids) - length(cell_id_list), cell_type);
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

edges = 0:bin_width:stim_dur;
n_bins = length(edges) - 1;

%% figure out what condition info is actually available in this dataset
combos = datarun.stimulus.combinations;
num_stim = length(combos);

has_direction = isfield(combos, 'DIRECTION');
has_tp        = isfield(combos, 'TEMPORAL_PERIOD');

if ~isempty(direction) && ~has_direction
    warning('saverasters:noDirectionField', ...
        '''direction'' filter given but this stimulus set has no DIRECTION field - ignoring filter.');
    direction = [];
end
if ~isempty(temporal_period) && ~has_tp
    warning('saverasters:noTPField', ...
        '''temporal_period'' filter given but this stimulus set has no TEMPORAL_PERIOD field - ignoring filter.');
    temporal_period = [];
end

%% determine which stimulus combinations to include
match = true(1, num_stim);
if has_direction && ~isempty(direction)
    all_dirs = [combos.DIRECTION];
    match = match & (all_dirs == direction);
end
if has_tp && ~isempty(temporal_period)
    all_tps = [combos.TEMPORAL_PERIOD];
    match = match & (all_tps == temporal_period);
end
stim_indices = find(match);

if isempty(stim_indices)
    error('saverasters:noMatch', 'No stimulus combinations matched the given filters.');
end

%% loop over cells, consolidating matched conditions into one matrix each
for c = 1:length(cell_id_list)
    cell_id  = cell_id_list(c);
    cell_idx = get_cell_indices(datarun, cell_id);
    spike_times = datarun.spikes{cell_idx};   % spike times (s) for this cell

    cell_name = matlab.lang.makeValidName(sprintf('cell%d', cell_id));

    all_blocks = cell(length(stim_indices), 1);
    col_names  = {};

    for k = 1:length(stim_indices)
        stim_index = stim_indices(k);
        this_combo = combos(stim_index);

        trial_inds = find(datarun.stimulus.trial_list == stim_index);
        trial_start_times = datarun.stimulus.triggers(trial_inds);
        n_trials = length(trial_inds);

        trial_matrix = zeros(n_bins, n_trials);
        for rep = 1:n_trials
            t0 = trial_start_times(rep);
            t1 = t0 + stim_dur;
            spk = spike_times(spike_times >= t0 & spike_times < t1) - t0;
            trial_matrix(:, rep) = histcounts(spk, edges)';
        end
        all_blocks{k} = trial_matrix;

        % build the condition label from whatever fields are actually available
        if has_direction && has_tp
            cond_name = sprintf('dir%d_TP%d', round(this_combo.DIRECTION), round(this_combo.TEMPORAL_PERIOD));
        elseif has_direction
            cond_name = sprintf('dir%d', round(this_combo.DIRECTION));
        elseif has_tp
            cond_name = sprintf('TP%d', round(this_combo.TEMPORAL_PERIOD));
        else
            cond_name = sprintf('cond%d', stim_index);
        end
        cond_name = matlab.lang.makeValidName(cond_name);

        for rep = 1:n_trials
            col_names{end+1} = sprintf('%s_trial%d', cond_name, rep); %#ok<AGROW>
        end
    end

    cell_matrix = horzcat(all_blocks{:});
    cell_table = array2table(cell_matrix, 'VariableNames', col_names);

    csv_path = fullfile(output_dir, [cell_name '.csv']);
    writetable(cell_table, csv_path);

    fprintf('Cell %d: saved consolidated matrix (%d bins x %d columns) to %s\n', ...
        cell_id, n_bins, size(cell_matrix, 2), csv_path);
end

fprintf('\nDone. %d cells saved to: %s\n', length(cell_id_list), output_dir);

end