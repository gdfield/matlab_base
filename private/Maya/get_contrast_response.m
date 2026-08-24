function crf = get_contrast_response(results, cell_ids, varying_field, varargin)
% GET_CONTRAST_RESPONSE   Compute a contrast-response function (mean firing
% rate vs. contrast) from condition-sorted spike times
%
% usage: crf = get_contrast_response(results, cell_ids, varying_field, varargin)
%
% arguments:  results        - struct array from GET_ALL_CONDITION_SPIKE_TIMES /
%                                  GET_ALL_CONDITION_SPIKE_TIMES1, one entry per
%                                  stimulus combination (.params, .spike_times)
%             cell_ids       - vector of cell ids to compute this for - must be
%                                  the SAME cell_ids, in the SAME order, that you
%                                  passed into get_all_condition_spike_times*
%                                  (results.spike_times is indexed by that order,
%                                  not by cell_id)
%             varying_field  - name of the field in results(k).params that
%                                  holds contrast, e.g. 'CONTRAST'. Field names
%                                  come straight from your stimulus file, and
%                                  this codebase doesn't standardize them - if
%                                  you're not sure of the exact name/casing,
%                                  run fieldnames(results(1).params) first.
%             varargin       - name-value pairs to fix every OTHER stimulus
%                                  parameter you want held constant, e.g.
%                                  'TEMPORAL_PERIOD', 60, 'SPATIAL_PERIOD', 160
%                                  Only combinations matching every fixed value
%                                  are included. Any parameter you don't fix
%                                  (e.g. direction) gets pooled together at each
%                                  contrast level rather than excluded - so fix
%                                  everything you don't want averaged over.
%
%                                  'stim_dur' is required (seconds) - use the
%                                  same value you passed into
%                                  get_all_condition_spike_times*, so firing
%                                  rate comes out in Hz rather than spikes/trial.
%
% outputs: crf - struct array, one entry per cell_id, in the same order as
% cell_ids:
%              .cell_id    - the cell id
%              .contrast   - sorted vector of contrast levels found
%              .mean_rate  - mean firing rate (Hz) at each contrast level
%              .sem_rate   - SEM across trials at each contrast level
%              .n_trials   - number of trials pooled at each contrast level
%
% example:
%
%   [results, datarun1] = get_all_condition_spike_times1(datarun1, mapped_cell_ids, 10, 'stim_dur', 8);
%   fieldnames(results(1).params)   % check the real field names once
%   crf = get_contrast_response(results, mapped_cell_ids, 'CONTRAST', ...
%       'stim_dur', 8, 'TEMPORAL_PERIOD', 60, 'SPATIAL_PERIOD', 160);
%   plot_contrast_response(crf)
%
% See also: GET_ALL_CONDITION_SPIKE_TIMES1, PLOT_CONTRAST_RESPONSE, BROWSE_CONTRAST_RESPONSE
%
% 2026-07 mcarleton
%

p = inputParser;
p.addParamValue('stim_dur', []);
p.KeepUnmatched = true;   % anything else is a fixed parameter to filter combos on
p.parse(varargin{:});
opts = p.Results;
fixed_params = p.Unmatched;

if isempty(opts.stim_dur)
    error(['get_contrast_response: ''stim_dur'' is required (seconds) - use the same value ' ...
        'you passed to get_all_condition_spike_times*, so firing rates come out in Hz.']);
end

if isempty(results)
    error('get_contrast_response: results is empty.');
end

fixed_names = fieldnames(fixed_params);

% FILTER: keep only combos that have varying_field and match every fixed param
matching = [];
for k = 1:length(results)
    combo_params = results(k).params;
    if ~isfield(combo_params, varying_field)
        continue
    end
    ok = true;
    for ff = 1:length(fixed_names)
        fname = fixed_names{ff};
        if ~isfield(combo_params, fname) || ...
                ~values_match(combo_params.(fname), fixed_params.(fname))
            ok = false;
            break
        end
    end
    if ok
        matching(end+1) = k; %#ok<AGROW>
    end
end

if isempty(matching)
    error(['get_contrast_response: no stimulus combinations matched. Check varying_field ' ...
        '(''%s'') and your fixed parameter names/values against fieldnames(results(1).params): %s'], ...
        varying_field, strjoin(fieldnames(results(1).params), ', '));
end

% collect the distinct contrast levels among matching combos
contrast_vals = nan(1, length(matching));
for mm = 1:length(matching)
    contrast_vals(mm) = to_numeric(results(matching(mm)).params.(varying_field));
end
unique_contrasts = unique(contrast_vals);

n_cells = length(cell_ids);
crf = struct('cell_id', cell(1,n_cells), 'contrast', cell(1,n_cells), ...
    'mean_rate', cell(1,n_cells), 'sem_rate', cell(1,n_cells), 'n_trials', cell(1,n_cells));

for c = 1:n_cells
    mean_rate = nan(size(unique_contrasts));
    sem_rate = nan(size(unique_contrasts));
    n_trials = zeros(size(unique_contrasts));

    for ci = 1:length(unique_contrasts)
        rates = [];
        for mm = 1:length(matching)
            k = matching(mm);
            if abs(to_numeric(results(k).params.(varying_field)) - unique_contrasts(ci)) > 1e-9
                continue
            end
            trial_spikes = results(k).spike_times(c,:);
            trial_rates = cellfun(@length, trial_spikes) / opts.stim_dur;
            rates = [rates, trial_rates]; %#ok<AGROW>
        end
        mean_rate(ci) = mean(rates);
        sem_rate(ci) = std(rates) / sqrt(length(rates));
        n_trials(ci) = length(rates);
    end

    crf(c).cell_id = cell_ids(c);
    crf(c).contrast = unique_contrasts;
    crf(c).mean_rate = mean_rate;
    crf(c).sem_rate = sem_rate;
    crf(c).n_trials = n_trials;
end

end


function v = to_numeric(x)
% coerce a param value (double, char, or string) to a double for comparison
if ischar(x) || isstring(x)
    v = str2double(x);
else
    v = double(x);
end
end


function tf = values_match(a, b)
% compare two param values loosely - numeric with tolerance if both parse as
% numbers, exact match otherwise (covers string-vs-numeric field storage)
an = to_numeric(a);
bn = to_numeric(b);
if ~isnan(an) && ~isnan(bn)
    tf = abs(an - bn) < 1e-9;
else
    tf = isequal(a, b);
end
end
