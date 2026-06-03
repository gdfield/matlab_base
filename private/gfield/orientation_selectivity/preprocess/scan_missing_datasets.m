%% scan_missing_datasets.m
%
% Scans the raw-data directory for the 23 dataset dates that currently
% have no _os.mat file in os_analysis/. For each dataset it lists:
%   - data<NNN> subfolders (both top-level and under Yass/)
%   - stimuli/ folder contents
%   - A best-guess pairing of grating + white-noise runs
%
% Output is printed to the command window AND written to a log file
% in the orientation_selectivity folder. Send that log file back to
% Claude and it will produce ready-to-paste os_datasets.m entries.
%
% Set `base_path` below to match where the raw data lives on your system.
%
% 2026-04 gdf / claude: Created

% ---- USER SETTINGS --------------------------------------------------
base_path = '/Volumes/gdf/rat-data/';
log_path  = fullfile(fileparts(mfilename('fullpath')), '..', ...
                     'missing_dataset_scan.txt');
% ---------------------------------------------------------------------

missing = {
    '2018-07-18-0'
    '2020-03-19-0'
    '2021-07-13-0'
    '2021-07-13-1'
    '2021-07-15-0'
    '2021-12-28-0'
    '2023-07-19-0'
    '2023-07-27-0'
    '2023-08-09-0'
    '2023-08-24-0'
    '2023-10-04-0'
    '2023-10-17-0'
    '2023-10-20-0'
    '2023-10-27-0'
    '2023-10-31-0'
    '2023-11-08-0'
    '2024-01-17-0'
    '2024-01-30-0'
    '2024-02-19-0'
    '2024-05-08-0'
    '2024-08-06-0'
    '2024-08-13-0'
    '2024-11-14-0'
};

fid = fopen(log_path, 'w');
assert(fid ~= -1, 'Could not open log file for writing: %s', log_path);

cleanup = onCleanup(@() fclose(fid));

log_both(fid, '==================================================\n');
log_both(fid, 'Missing-dataset scan\n');
log_both(fid, 'base_path = %s\n', base_path);
log_both(fid, 'date      = %s\n', datestr(now));
log_both(fid, '==================================================\n\n');

for k = 1:numel(missing)
    ds = missing{k};
    ds_root = fullfile(base_path, ds);

    log_both(fid, '---- %s ----\n', ds);
    log_both(fid, '  path: %s\n', ds_root);

    if ~isfolder(ds_root)
        log_both(fid, '  [MISSING] folder not found on disk\n\n');
        continue
    end

    % 1. list data<NNN> folders at top level and under Yass/
    top_data  = list_data_dirs(ds_root);
    yass_root = fullfile(ds_root, 'Yass');
    if isfolder(yass_root)
        yass_data = list_data_dirs(yass_root);
    else
        yass_data = {};
    end

    log_both(fid, '  top-level data folders (%d):\n', numel(top_data));
    for i = 1:numel(top_data)
        log_both(fid, '    %s\n', top_data{i});
    end
    if ~isempty(yass_data)
        log_both(fid, '  Yass/ data folders (%d):\n', numel(yass_data));
        for i = 1:numel(yass_data)
            log_both(fid, '    Yass/%s\n', yass_data{i});
        end
    end

    % 2. list stimuli/ folder contents
    stim_root = fullfile(ds_root, 'stimuli');
    if isfolder(stim_root)
        d = dir(stim_root);
        d = d(~[d.isdir]);
        log_both(fid, '  stimuli/ files (%d):\n', numel(d));
        for i = 1:numel(d)
            log_both(fid, '    %s\n', d(i).name);
        end
    else
        log_both(fid, '  [no stimuli/ folder]\n');
    end

    % 3. heuristic pairing suggestion
    %    A grating run is a dataNNN that has a matching stimulus file
    %    (sNN.txt, dataNNN.txt, dgNNN.txt, etc.). WN run is typically an
    %    adjacent dataNNN (usually N-1) with no stimulus file.
    log_both(fid, '  --- suggested pairings ---\n');
    all_data = [strcat('top/', top_data); strcat('Yass/', yass_data)];
    suggest_pairs(fid, stim_root, all_data);

    log_both(fid, '\n');
end

log_both(fid, 'Log written to:\n  %s\n', log_path);
fprintf('\nDone. Please send the log file back:\n  %s\n', log_path);

% =====================================================================
% helpers
% =====================================================================
function log_both(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

function names = list_data_dirs(root)
    d = dir(root);
    d = d([d.isdir]);
    names = {};
    for i = 1:numel(d)
        n = d(i).name;
        if ~isempty(regexp(n, '^data\d+', 'once'))
            names{end+1,1} = n; %#ok<AGROW>
        end
    end
    names = sort(names);
end

function suggest_pairs(fid, stim_root, all_data)
    % pull the leading "dataNNN" number out of each data folder name
    nums = nan(numel(all_data),1);
    for i = 1:numel(all_data)
        folder_only = all_data{i};
        slash = strfind(folder_only, '/');
        folder_only = folder_only(slash(end)+1:end);
        t = regexp(folder_only, '^data(\d+)', 'tokens', 'once');
        if ~isempty(t), nums(i) = str2double(t{1}); end
    end

    % list of stimulus file basenames
    stim_files = {};
    if isfolder(stim_root)
        d = dir(stim_root);
        stim_files = {d(~[d.isdir]).name};
    end

    % flag dataNNN as "grating-like" if a matching stim file exists
    is_grating = false(numel(all_data),1);
    stim_match = cell(numel(all_data),1);
    for i = 1:numel(all_data)
        if isnan(nums(i)), continue; end
        nnn = sprintf('%03d', nums(i));
        nn  = sprintf('%02d',  mod(nums(i), 100));  % e.g. 07 for data007
        candidates = {['data', nnn, '.txt'], ...
                      ['s',     nn,  '.txt'], ['s',  nn],  ...
                      ['dg',   nnn,  '.txt'], ['dg', nn, '.txt'], ...
                      ['S',     nn], ['S', nn, '.txt']};
        for c = 1:numel(candidates)
            if any(strcmpi(stim_files, candidates{c}))
                is_grating(i) = true;
                stim_match{i} = candidates{c};
                break
            end
        end
    end

    % for each grating candidate, suggest nearest-lower dataNNN as WN
    if ~any(is_grating)
        fprintf(fid, '    (no dataNNN folder matched a stimuli/ file -- manual review needed)\n');
        fprintf('    (no dataNNN folder matched a stimuli/ file -- manual review needed)\n');
        return
    end

    for i = 1:numel(all_data)
        if ~is_grating(i), continue; end
        % nearest-lower number among remaining
        diffs = nums(i) - nums;
        diffs(diffs <= 0 | isnan(diffs)) = inf;
        [~, j] = min(diffs);
        if isinf(diffs(j))
            wn_guess = '(none found)';
        else
            wn_guess = all_data{j};
        end
        fprintf(fid, '    grating: %s   (stim=%s)   ->   WN guess: %s\n', ...
                all_data{i}, stim_match{i}, wn_guess);
        fprintf('    grating: %s   (stim=%s)   ->   WN guess: %s\n', ...
                all_data{i}, stim_match{i}, wn_guess);
    end
end
