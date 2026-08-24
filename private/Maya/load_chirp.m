function datarun = load_chirp(datarun, chirp_path, varargin)
% LOAD_CHIRP   Load precomputed chirp-stimulus analysis output into datarun
%
% usage: datarun = load_chirp(datarun, chirp_path, varargin)
%
% arguments:  datarun    - datarun struct. datarun.cell_ids must already exist
%                             (i.e. run after load_neurons / load_params)
%             chirp_path - full path to the {algorithm}_{search}.mat file written
%                             by the chirp analysis pipeline, e.g.
%                             '.../kilosort4/kilosort4_ChirpStimulus.mat'
%             varargin   - struct or list of optional parameters (see below)
%
% outputs: datarun with a new datarun.chirp field, indexed by cell_number the
% same way datarun.ei.eis / datarun.stas.stas are (i.e. datarun.chirp.psth_mean{
% get_cell_indices(datarun,cell_id) }), so it plugs into the rest of this
% codebase's conventions:
%
%   datarun.chirp.psth_mean{cn}      - (n_bins,) Hz, trial-averaged PSTH
%   datarun.chirp.spikes_binned{cn}  - (n_trials, n_bins) raw per-trial spike counts
%   datarun.chirp.quality_index(cn)  - scalar QI, NaN if the cell was silent
%   datarun.chirp.responsive(cn)     - logical, QI > qi_threshold (NaNs -> false)
%   datarun.chirp.matched(cn)        - logical, true if this cell_id had a row in
%                                         the chirp file's cluster_id list
%   datarun.chirp.phase.step_on/step_off/freq_sweep/contrast
%                                     - [start end] bin indices, 1-indexed and
%                                         end-INCLUSIVE (MATLAB slicing convention -
%                                         see note #2 below for how this differs
%                                         from the reference doc's Python examples)
%   datarun.chirp.params             - struct: bin_size_ms, pre_ms, stim_ms, step_ms,
%                                         inter_ms, freq_ms, contrast_ms, freq_min_hz,
%                                         freq_max_hz
%   datarun.chirp.source             - struct: experiment_name, chunk_name,
%                                         sort_algorithm, source_files, file
%
% cells in datarun.cell_ids with no match in the chirp file are left
% empty/NaN/false, so every array here is always length(datarun.cell_ids) long -
% safe to index with the usual cell_number from get_cell_indices.
%
% optional params, their default values, and what they specify:
%
% qi_threshold      1       QI above this counts as 'responsive' (matches the
%                               rule of thumb in the chirp reference doc)
% verbose           true    print a match/responsiveness summary after loading
%
%
% ***** VERIFY THESE BEFORE TRUSTING RESULTS *****
%
% 1. ID CONVENTION - this assumes the chirp file's cluster_id values are the SAME
%    IDs as datarun.cell_ids (e.g. because the chirp stimulus was sorted in the
%    same kilosort chunk as this datarun). If the chirp analysis was run as a
%    separate sort/chunk, cluster_id will NOT line up with datarun.cell_ids and
%    the match will silently come back mostly/entirely empty. Check the printed
%    "matched N / M cells" summary the first time you run this - if N is much
%    smaller than M, you likely need to load the chirp stimulus's own data run
%    separately and use MAP_EI to relate its cells to this datarun instead
%    (see the gratings section of MCpipeline.m for the pattern).
%
% 2. BIN INDEXING - phase_* fields in the .mat file were computed in Python with
%    0-indexed, end-exclusive slicing (psth_mean[:, start:end] means bins
%    start..end-1). This function converts to 1-indexed, end-inclusive MATLAB
%    ranges, i.e. datarun.chirp.psth_mean{cn}(s:e) with
%    [s, e] = datarun.chirp.phase.freq_sweep should give you the identical set
%    of bins as the Python start:end example in the reference doc. Sanity-check
%    this once against a real file (e.g. phase.step_on(1) should correspond to
%    roughly pre_ms / bin_size_ms).
%
% 2026-07 mcarleton
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addParamValue('qi_threshold', 1);
p.addParamValue('verbose', true);
p.parse(varargin{:});
params = p.Results;


% CHECK INPUTS

if ~isfield(datarun,'cell_ids') || isempty(datarun.cell_ids)
    error('load_chirp: datarun.cell_ids does not exist yet - run load_neurons/load_params first.');
end

if ~exist('chirp_path','var') || isempty(chirp_path)
    error('load_chirp: chirp_path is required, e.g. .../kilosort4_ChirpStimulus.mat');
end

if ~exist(chirp_path,'file')
    error('load_chirp: no file found at ''%s''', chirp_path);
end


% LOAD THE FILE

d = load(chirp_path);

required_fields = {'cluster_id','psth_mean','spikes_binned','quality_index'};
for rr = 1:length(required_fields)
    if ~isfield(d, required_fields{rr})
        error('load_chirp: expected field ''%s'' not found in %s', required_fields{rr}, chirp_path);
    end
end

chirp_cluster_ids = double(d.cluster_id(:));
n_cells_total = length(datarun.cell_ids);


% MATCH CHIRP ROWS TO datarun.cell_ids, IN cell_number ORDER

datarun.chirp.psth_mean = cell(n_cells_total, 1);
datarun.chirp.spikes_binned = cell(n_cells_total, 1);
datarun.chirp.quality_index = nan(n_cells_total, 1);
datarun.chirp.matched = false(n_cells_total, 1);

n_matched = 0;
for cn = 1:n_cells_total
    row = find(chirp_cluster_ids == datarun.cell_ids(cn), 1);
    if isempty(row)
        continue
    end
    datarun.chirp.psth_mean{cn} = d.psth_mean(row,:);
    datarun.chirp.spikes_binned{cn} = squeeze(d.spikes_binned(row,:,:));
    datarun.chirp.quality_index(cn) = d.quality_index(row);
    datarun.chirp.matched(cn) = true;
    n_matched = n_matched + 1;
end

% NaN = silent cell (denominator was 0 when QI was computed) - filter with
% isnan, not with > 0, per the reference doc
datarun.chirp.responsive = datarun.chirp.matched & ~isnan(datarun.chirp.quality_index) ...
    & (datarun.chirp.quality_index > params.qi_threshold);


% TIMING PARAMS (SCALAR, STIMULUS-WIDE - SAME FOR EVERY CELL)

timing_fields = {'bin_size_ms','pre_ms','stim_ms','step_ms','inter_ms','freq_ms', ...
    'contrast_ms','freq_min_hz','freq_max_hz'};
datarun.chirp.params = struct;
for tt = 1:length(timing_fields)
    if isfield(d, timing_fields{tt})
        val = d.(timing_fields{tt});
        if numel(val) == 1
            val = double(val);
        end
        datarun.chirp.params.(timing_fields{tt}) = val;
    end
end


% PHASE BIN INDICES - Python 0-indexed/end-exclusive -> MATLAB 1-indexed/end-inclusive
% Python  [start, end)  ==  MATLAB  (start+1):end

phase_fields = {'phase_step_on','phase_step_off','phase_freq_sweep','phase_contrast'};
phase_names  = {'step_on','step_off','freq_sweep','contrast'};
datarun.chirp.phase = struct;
for pp = 1:length(phase_fields)
    if isfield(d, phase_fields{pp})
        se = double(d.(phase_fields{pp})(:))';
        datarun.chirp.phase.(phase_names{pp}) = [se(1)+1, se(2)];
    end
end


% BOOKKEEPING / PROVENANCE

source_fields = {'experiment_name','chunk_name','sort_algorithm','source_files'};
datarun.chirp.source = struct;
for ss = 1:length(source_fields)
    if isfield(d, source_fields{ss})
        datarun.chirp.source.(source_fields{ss}) = d.(source_fields{ss});
    end
end
datarun.chirp.source.file = chirp_path;


% SUMMARY

if params.verbose
    fprintf('\nload_chirp: matched %d / %d cells to %s\n', n_matched, n_cells_total, chirp_path);
    if n_matched > 0
        fprintf('load_chirp: %d / %d matched cells pass QI > %g (responsive)\n', ...
            sum(datarun.chirp.responsive), n_matched, params.qi_threshold);
    end
    if n_matched < n_cells_total
        fprintf(['load_chirp: %d cell_ids had no match in cluster_id. If this number is large,\n' ...
            '   the chirp file was probably sorted separately from this datarun - see note #1\n' ...
            '   in this function''s header before trusting datarun.chirp.\n'], n_cells_total - n_matched);
    end
end
