%% os_cell_finder_script_claude_auto.m
%
% Non-interactive batch version of os_cell_finder_script.m for use
% via automated MATLAB sessions (no stdin available).
%
% What this script does vs. the interactive version:
%   - Processes the datasets in new_os_datasets.m (optionally a subset via batch_override)
%   - Skips white-noise RF mapping entirely
%   - Exports one multi-panel PDF per OS-candidate cell (rasters + polar)
%   - Auto-saves all threshold-passing cells to os_maybe_list; os_cell_list
%     is left empty for manual review by the experimenter after inspecting PDFs
%
% THRESHOLD CONSTANTS -- keep in sync with os_cell_finder_script.m:
%   OSI_thresh, DSI_thresh, response_cor_thresh
%
% Output locations:
%   PDFs:  <os_analysis_root>/claude_auto_plots/<dataset>/cell<id>_sp<sp>_tp<tp>.pdf
%   .mat:  <os_analysis_root>/<dataset>_os.mat
%
% Usage: run this script from the MATLAB command window or via script runner.
%   Results are printed to the command window and written to disk.
%
% 2026-09-26 generated for OS-RGC manuscript curation pipeline

%% --- Parameters (keep in sync with os_cell_finder_script.m) ---
OSI_thresh          = 0.30;
DSI_thresh          = 0.25;
response_cor_thresh = 0.30;

%% --- Paths ---
os_analysis_root = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
plot_root        = [os_analysis_root, 'claude_auto_plots/'];

% Suppress figure windows during batch run
set(0, 'DefaultFigureVisible', 'off');

%% --- Load dataset list ---
data_list = new_os_datasets();
if exist('batch_override', 'var') && ~isempty(batch_override)   % optional dataset subset
    keep = cellfun(@(p) any(contains(p, batch_override)), {data_list.grating_datapath});
    data_list = data_list(keep);
end
num_datasets = length(data_list);

%% --- Summary accumulators ---
summary = struct();

fprintf('\n========================================\n');
fprintf('OS-RGC batch curation: %d datasets\n', num_datasets);
fprintf('OSI_thresh=%.2f  DSI_thresh=%.2f  corr_thresh=%.2f\n', ...
    OSI_thresh, DSI_thresh, response_cor_thresh);
fprintf('========================================\n\n');

%% --- Main loop ---
for dset = 1:num_datasets

    dl = data_list(dset);

    % Extract dataset name from grating path (4th path component)
    split_path  = regexp(dl.grating_datapath, '/', 'split');
    dataset_name = split_path{5};  % e.g. '2023-10-04-0'

    fprintf('[%d/%d] %s\n', dset, num_datasets, dataset_name);

    %% Load grating datarun
    try
        datarun = load_data(dl.grating_datapath);
        datarun = load_neurons(datarun);
        datarun = load_params(datarun);
    catch ME
        fprintf('  ERROR loading datarun: %s\n', ME.message);
        summary(dset).dataset     = dataset_name;
        summary(dset).status      = 'load_error';
        summary(dset).n_candidate = 0;
        summary(dset).cell_ids    = [];
        continue
    end

    %% Load stimulus via trigger check
    try
        datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
    catch ME
        fprintf('  ERROR loading stimulus: %s\n', ME.message);
        summary(dset).dataset     = dataset_name;
        summary(dset).status      = 'stim_error';
        summary(dset).n_candidate = 0;
        summary(dset).cell_ids    = [];
        continue
    end

    % Verify stimulus loaded with usable repetitions
    if ~isfield(datarun, 'stimulus') || isempty(datarun.stimulus)
        fprintf('  SKIP: stimulus did not load\n');
        summary(dset).dataset     = dataset_name;
        summary(dset).status      = 'no_stimulus';
        summary(dset).n_candidate = 0;
        summary(dset).cell_ids    = [];
        continue
    end

    spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
    temp_periods    = datarun.stimulus.params.TEMPORAL_PERIOD;
    num_sps = length(spatial_periods);
    num_tps = length(temp_periods);
    num_dirs = length(datarun.stimulus.params.DIRECTION);
    num_rgcs = length(datarun.cell_ids);

    fprintf('  %d RGCs  |  %d dirs  |  %d SPs  |  %d TPs\n', ...
        num_rgcs, num_dirs, num_sps, num_tps);

    %% Compute OSI / DSI / reliability for all cells
    try
        [OSIs, DSIs, mean_corr] = find_OS_RGCs(datarun, 'all', ...
            'plot_tuning', false, 'verbose', false, 'pause', false);
    catch ME
        fprintf('  ERROR in find_OS_RGCs: %s\n', ME.message);
        summary(dset).dataset     = dataset_name;
        summary(dset).status      = 'osi_error';
        summary(dset).n_candidate = 0;
        summary(dset).cell_ids    = [];
        continue
    end

    %% Create output directories
    dset_plot_dir = [plot_root, dataset_name, '/'];
    if ~exist(dset_plot_dir, 'dir')
        mkdir(dset_plot_dir);
    end

    %% Identify OS candidates (computed and saved before any plotting)
    os_cell_list  = [];   % left empty -- experimenter promotes from maybe after PDF review
    os_maybe_list = [];
    cand_rows     = [];
    for rgc = 1:num_rgcs
        insig_DSI = all(DSIs(rgc,:,:)   < DSI_thresh,       'all');
        sig_OSIs  = any(OSIs(rgc,:,:)   > OSI_thresh,        'all');
        sig_cors  = any(mean_corr(rgc,:,:) > response_cor_thresh, 'all');
        if sig_OSIs && insig_DSI && sig_cors
            os_maybe_list = [os_maybe_list, datarun.cell_ids(rgc)]; %#ok<AGROW>
            cand_rows     = [cand_rows, rgc]; %#ok<AGROW>
        end
    end
    n_maybe = length(os_maybe_list);

    %% Save _os.mat
    mat_path = [os_analysis_root, dataset_name, '_os.mat'];
    save(mat_path, 'os_cell_list', 'os_maybe_list');
    fprintf('  OS candidates: %d / %d  -->  saved %s\n', n_maybe, num_rgcs, mat_path);

    %% Export plots: one PDF per candidate per SP x TP combination
    % optional workspace flags: batch_skip_existing_pdfs (resume), batch_plot_budget_s
    skip_existing = exist('batch_skip_existing_pdfs', 'var') && batch_skip_existing_pdfs;
    if exist('batch_plot_budget_s', 'var'), budget = batch_plot_budget_s; else, budget = Inf; end
    t_plot = tic; plots_complete = true;
    fig_handle = figure('Visible', 'off');
    for rgc = cand_rows
        cell_id = datarun.cell_ids(rgc);
        for spat_p = 1:num_sps
            for temp_p = 1:num_tps
                temp_sp = spatial_periods(spat_p);
                temp_tp = temp_periods(temp_p);
                pdf_name = sprintf('cell%d_sp%g_tp%g.pdf', cell_id, temp_sp, temp_tp);
                if skip_existing && exist([dset_plot_dir, pdf_name], 'file'), continue; end
                if toc(t_plot) > budget, plots_complete = false; break; end
                fig_title_str = sprintf('%s  cell %d  SP %.0f  TP %.0f  OSI %.2f', ...
                    dataset_name, cell_id, temp_sp, temp_tp, max(OSIs(rgc, spat_p, temp_p)));
                try
                    [temp_tuning, temp_spike_nums, ~] = get_direction_tuning( ...
                        datarun.spikes{rgc}, datarun.stimulus, 'SP', temp_sp, 'TP', temp_tp);
                    plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, ...
                        'fig_num', fig_handle, 'fig_title', fig_title_str, 'clear_fig', true);
                    drawnow
                    exportgraphics(fig_handle, [dset_plot_dir, pdf_name], 'ContentType', 'image');
                catch ME
                    fprintf('  WARNING: plot failed for cell %d SP %g TP %g: %s\n', ...
                        cell_id, temp_sp, temp_tp, ME.message);
                end
            end
        end
    end
    close(fig_handle);
    if ~plots_complete
        fprintf('  PLOTS INCOMPLETE (time budget): rerun with batch_skip_existing_pdfs = true\n');
    end

    summary(dset).dataset     = dataset_name;
    summary(dset).status      = 'ok';
    summary(dset).n_candidate = n_maybe;
    summary(dset).cell_ids    = os_maybe_list;

end  % dataset loop

%% Re-enable figure visibility
set(0, 'DefaultFigureVisible', 'on');

%% Print final summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n');
fprintf('%-20s  %-12s  %s\n', 'Dataset', 'Status', 'Candidates');
for d = 1:length(summary)
    if ~isempty(summary(d).dataset)
        fprintf('%-20s  %-12s  %d\n', ...
            summary(d).dataset, summary(d).status, summary(d).n_candidate);
    end
end
fprintf('========================================\n');
