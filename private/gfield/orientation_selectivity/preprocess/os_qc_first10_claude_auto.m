%% os_qc_first10_claude_auto.m
%
% QC scan: plots rasters + polar tuning for the first 10 cells (by index,
% NOT by OSI) in each of the low-candidate-count datasets from the
% 2026-09-26/27 batch curation run. Purpose: sanity-check that the
% grating stimulus is being parsed correctly (reasonable, time-locked
% responses) independent of whether cells clear the OS threshold.
%
% Datasets included here all returned fewer than 20 OS candidates in the
% batch run:
%   2021-07-13-0   3 candidates / 322 cells
%   2021-07-13-1   3 candidates / 285 cells
%   2021-07-15-0   9 candidates / 369 cells
%   2021-12-28-0  11 candidates / 362 cells
%   2023-10-27-0   9 candidates / 400 cells
%   2024-01-17-0   0 candidates / 483 cells
%   2024-08-06-0   3 candidates / 487 cells  (load_stim trimmed 6 reps -> 5)
%   2024-08-13-0   0 candidates / 506 cells  (load_stim trimmed 6 reps -> 5)
%   2024-11-14-0   6 candidates / 690 cells  (load_stim trigger count off by 1)
%
% The last three above also printed a load_stim/load_stim_amr trigger-count
% mismatch warning during the batch run, so they deserve an extra careful
% look, that warning could mean a dropped trial rather than a real absence
% of OS cells.
%
% By default this plots ONE representative SP/TP combination per cell
% (the first entry in each), to keep output at 10 PDFs per dataset. Set
% plot_all_combos = true below to export every SP x TP combination
% instead (more PDFs, more complete picture).
%
% Output: <os_analysis_root>/claude_auto_plots/qc_first10/<dataset>/cell<id>_sp<sp>_tp<tp>.pdf
%
% 2026-09-27 generated for OS-RGC manuscript curation pipeline QC pass

%% --- Parameters ---
n_cells_to_check = 10;
plot_all_combos  = false;   %#ok<NASGU> % set true to export every SP x TP combo per cell

qc_dataset_names = {
    '2021-07-13-0', '2021-07-13-1', '2021-07-15-0', '2021-12-28-0', ...
    '2023-10-27-0', '2024-01-17-0', '2024-08-06-0', '2024-08-13-0', ...
    '2024-11-14-0'};

% optional: set qc_override = {'<dataset>', ...} in the workspace before
% running to QC only those datasets (clear it afterward to restore the full list)
if exist('qc_override', 'var') && ~isempty(qc_override)
    qc_dataset_names = qc_override;
end

%% --- Paths ---
os_analysis_root = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
plot_root        = [os_analysis_root, 'claude_auto_plots/qc_first10/'];

set(0, 'DefaultFigureVisible', 'off');

%% --- Load full dataset registry and filter to the QC subset ---
all_datasets = new_os_datasets();

qc_list = [];
for k = 1:numel(all_datasets)
    split_path = regexp(all_datasets(k).grating_datapath, '/', 'split');
    dname = split_path{5};
    if ismember(dname, qc_dataset_names)
        if isempty(qc_list)
            qc_list = all_datasets(k);
        else
            qc_list(end+1) = all_datasets(k); %#ok<AGROW>
        end
    end
end
fprintf('QC scan: %d of %d requested datasets found in registry\n', numel(qc_list), numel(qc_dataset_names));

%% --- Main loop ---
fig_handle = figure('Visible', 'off');

for dset = 1:numel(qc_list)

    dl = qc_list(dset);
    split_path   = regexp(dl.grating_datapath, '/', 'split');
    dataset_name = split_path{5};

    fprintf('\n[%d/%d] %s\n', dset, numel(qc_list), dataset_name);

    try
        datarun = load_data(dl.grating_datapath);
        datarun = load_neurons(datarun);
        datarun = load_params(datarun);
    catch ME
        fprintf('  ERROR loading datarun: %s\n', ME.message);
        continue
    end

    try
        datarun = os_trigger_check(datarun, dl.stimulus_path, dl.trigger_interval);
    catch ME
        fprintf('  ERROR loading stimulus: %s\n', ME.message);
        continue
    end

    if ~isfield(datarun, 'stimulus') || isempty(datarun.stimulus)
        fprintf('  SKIP: stimulus did not load\n');
        continue
    end

    spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
    temp_periods    = datarun.stimulus.params.TEMPORAL_PERIOD;
    num_sps = length(spatial_periods);
    num_tps = length(temp_periods);
    num_dirs = length(datarun.stimulus.params.DIRECTION);
    num_rgcs = length(datarun.cell_ids);

    fprintf('  %d RGCs total | %d dirs | %d SPs | %d TPs | checking first %d cells\n', ...
        num_rgcs, num_dirs, num_sps, num_tps, min(n_cells_to_check, num_rgcs));

    dset_plot_dir = [plot_root, dataset_name, '/'];
    if ~exist(dset_plot_dir, 'dir')
        mkdir(dset_plot_dir);
    end

    n_check = min(n_cells_to_check, num_rgcs);
    if plot_all_combos
        sp_range = 1:num_sps;
        tp_range = 1:num_tps;
    else
        sp_range = 1;
        tp_range = 1;
    end

    for rgc = 1:n_check

        cell_id = datarun.cell_ids(rgc);
        temp_spike_times = datarun.spikes{rgc};

        for spat_p = sp_range
            for temp_p = tp_range

                temp_sp = spatial_periods(spat_p);
                temp_tp = temp_periods(temp_p);

                fig_title_str = sprintf('%s  cell %d  (index %d/%d)  SP %.0f  TP %.0f', ...
                    dataset_name, cell_id, rgc, n_check, temp_sp, temp_tp);

                try
                    [temp_tuning, temp_spike_nums, ~] = get_direction_tuning( ...
                        temp_spike_times, datarun.stimulus, ...
                        'SP', temp_sp, 'TP', temp_tp);

                    plot_direction_tuning(temp_tuning, temp_spike_nums, ...
                        datarun.stimulus, ...
                        'fig_num',   fig_handle, ...
                        'fig_title', fig_title_str, ...
                        'clear_fig', true);
                    drawnow

                    pdf_name = sprintf('cell%d_sp%g_tp%g.pdf', cell_id, temp_sp, temp_tp);
                    exportgraphics(fig_handle, [dset_plot_dir, pdf_name], ...
                        'ContentType', 'image');

                catch ME
                    fprintf('  WARNING: plot failed for cell %d SP %g TP %g: %s\n', ...
                        cell_id, temp_sp, temp_tp, ME.message);
                end
            end
        end
    end

    fprintf('  Done: %s\n', dset_plot_dir);
end

close(fig_handle);
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nQC scan complete. PDFs in: %s\n', plot_root);
