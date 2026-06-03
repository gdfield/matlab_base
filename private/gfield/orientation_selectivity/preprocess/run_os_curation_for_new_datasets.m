%% run_os_curation_for_new_datasets.m
%
% Wrapper around the os_cell_finder_script curation workflow that runs
% only on a specified subset of dataset indices from os_datasets(), and
% writes the resulting _os.mat files into the os_analysis/ folder that
% already sits under private/gfield/orientation_selectivity/.
%
% Usage: edit `dset_indices` below to match the indices (in os_datasets.m)
% of the datasets you want to curate, then run the script. For each
% candidate cell you will be prompted for y (OS), n (not OS), or m (maybe).
%
% 2026-04 gdf / claude: Created

% ---- USER SETTINGS --------------------------------------------------
dset_indices = [];   % <-- set this, e.g. [32 35 37 38]

% location where _os.mat files should land (matches existing files)
os_out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'os_analysis');

% thresholds (kept identical to os_cell_finder_script.m)
OSI_thresh          = 0.30;
DSI_thresh          = 0.25;
response_cor_thresh = 0.30;
ei_corr_threshold   = 0.90;
filt_rad            = 0.75;
fast_plot           = true;
wn_map              = 1;
% ---------------------------------------------------------------------

assert(~isempty(dset_indices), ...
    'Set dset_indices at the top of this script before running.');

if ~isfolder(os_out_dir), mkdir(os_out_dir); end

data_list = os_datasets('experiment', 'all');

for kk = 1:numel(dset_indices)
    dset = dset_indices(kk);
    fprintf('\n=============================================\n');
    fprintf(' Curating dataset index %d  (%s)\n', dset, ...
            regexprep(data_list(dset).grating_datapath, '.*/(\d{4}-\d{2}-\d{2}-\d)/.*', '$1'));
    fprintf('=============================================\n');

    % --- load grating dataset
    datarun = load_data(data_list(dset).grating_datapath);
    datarun = load_neurons(datarun);
    datarun = load_params(datarun);
    datarun = load_ei(datarun, 'all');

    datarun = os_trigger_check(datarun, ...
                               data_list(dset).stimulus_path, ...
                               data_list(dset).trigger_interval);

    spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
    temp_periods    = datarun.stimulus.params.TEMPORAL_PERIOD;
    num_sps         = numel(spatial_periods);
    num_tps         = numel(temp_periods);

    % --- load WN dataset
    datarun_wn = load_data(data_list(dset).wn_datapath);
    datarun_wn = load_neurons(datarun_wn);
    datarun_wn = load_params(datarun_wn);
    datarun_wn = load_sta(datarun_wn, 'load_sta', 'all');
    datarun_wn = load_ei(datarun_wn, 'all');

    % --- run OS detection on all cells
    num_rgcs = numel(datarun.cell_ids);
    [OSIs, DSIs, mean_corr] = find_OS_RGCs(datarun, 'all', ...
                                   'plot_tuning', false, ...
                                   'verbose', false, 'pause', false);

    os_cell_list  = [];
    os_maybe_list = [];

    for rgc = 1:num_rgcs
        insig_DSI = all(DSIs(rgc,:,:) < DSI_thresh);
        sig_OSIs  = any(OSIs(rgc,:,:) > OSI_thresh);
        sig_cors  = any(mean_corr(rgc,:,:) > response_cor_thresh);

        if ~(sig_OSIs && insig_DSI && sig_cors), continue; end

        % --- plot direction-tuning rasters for each SP/TP
        fig_counter = 1;
        for sp = 1:num_sps
            for tp = 1:num_tps
                tspk = datarun.spikes{rgc};
                ftitle = sprintf('cell %d  SP %g  TP %g', ...
                                 datarun.cell_ids(rgc), ...
                                 spatial_periods(sp), temp_periods(tp));
                [tun, sn, ~] = get_direction_tuning(tspk, datarun.stimulus, ...
                                  'SP', spatial_periods(sp), ...
                                  'TP', temp_periods(tp));
                plot_direction_tuning(tun, sn, datarun.stimulus, ...
                                  'fig_num', fig_counter, ...
                                  'fig_title', ftitle, ...
                                  'fast_plot', fast_plot);
                fig_counter = fig_counter + 1;
            end
        end
        fprintf('OSIs %g  DSIs %g  corrs %g\n', ...
                OSIs(rgc,:,:), DSIs(rgc,:,:), mean_corr(rgc,:,:));

        % --- try to map to WN and plot RF
        if wn_map
            mapped = map_ei(datarun, datarun_wn, ...
                            'master_cell_type', datarun.cell_ids(rgc), ...
                            'corr_threshold',   ei_corr_threshold);
            idx = get_cell_indices(datarun, datarun.cell_ids(rgc));
            wn_id = mapped{idx};
            figure(10); clf;
            if isempty(wn_id)
                warning('RGC %d failed to map to WN', datarun.cell_ids(rgc));
            else
                rf = get_rf(datarun_wn, wn_id);
                [frf, ~] = rf_filtered(rf, struct('filt_type','gauss', ...
                                                  'radius',filt_rad));
                imagesc(squeeze(norm_image(frf(:,:,1))));
                colormap(brewermap([],'RdBu')); caxis([0,1]); axis equal;
                title(num2str(wn_id));
            end
        end

        % --- user prompt
        while true
            r = input('keep this cell? y/n/m:', 's');
            if strcmp(r,'y')
                os_cell_list(end+1) = datarun.cell_ids(rgc); %#ok<SAGROW>
                break
            elseif strcmp(r,'n')
                break
            elseif strcmp(r,'m')
                os_maybe_list(end+1) = datarun.cell_ids(rgc); %#ok<SAGROW>
                break
            else
                warning('enter y, n, or m');
            end
        end
    end

    % --- save
    split_path = regexp(data_list(dset).grating_datapath, '/', 'split');
    date_token = '';
    for s = 1:numel(split_path)
        if ~isempty(regexp(split_path{s}, '^\d{4}-\d{2}-\d{2}-\d$', 'once'))
            date_token = split_path{s}; break
        end
    end
    assert(~isempty(date_token), 'Could not extract date token from path');
    out_file = fullfile(os_out_dir, [date_token, '_os.mat']);
    save(out_file, 'os_cell_list', 'os_maybe_list');
    fprintf('Saved %s  (%d confident, %d maybe)\n', ...
            out_file, numel(os_cell_list), numel(os_maybe_list));
end

fprintf('\nAll requested datasets curated.\n');
