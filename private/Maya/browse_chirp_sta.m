function browse_chirp_sta(datarun, datarun1, master_ids_matched, mapped_cell_ids, varargin)
%BROWSE_CHIRP_STA  Step through matched cells, plotting the STA time
%course (from datarun) next to the chirp PSTH (from datarun1), advancing
%to the next pair when you click "Next".
%
%   browse_chirp_sta(datarun, datarun1, master_ids_matched, mapped_cell_ids)
%
% INPUTS
%   datarun             - master/WN datarun. Needs datarun.stas.time_courses
%                            populated - this is NOT done by load_sta itself,
%                            run this once first:
%                               datarun = get_sta_summaries(datarun, cell_spec);
%                            (cell_spec can just be master_ids_matched, or
%                            'all' / a cell type - see GET_STA_SUMMARIES)
%   datarun1            - slave/mapped datarun, with datarun1.chirp already
%                            populated (see LOAD_CHIRP) - chirp panel comes
%                            from here
%   master_ids_matched  - cell ids in datarun (time course side). In
%                            MCpipeline.m this is:
%                            datarun.cell_ids(master_idx_matched)
%   mapped_cell_ids     - cell ids in datarun1 (chirp side). In MCpipeline.m
%                            this is: cell2mat(cell_map)
%
%   master_ids_matched(k) and mapped_cell_ids(k) must already be paired 1:1
%   in matching order. This holds automatically for the master_ids_matched/
%   mapped_cell_ids built in MCpipeline.m, since both are derived by walking
%   cell_map (from MAP_EI) in the same left-to-right order - cell2mat just
%   drops the empty (unmatched) entries, and find(~cellfun(@isempty,...))
%   returns the surviving indices in that same order. If you build these
%   lists some other way, make sure that pairing still holds before using
%   this function - it does NOT re-match master to slave itself.
%
%   varargin - optional parameters (see below)
%
% optional params, their default values, and what they specify:
%
% figure        0       figure number to plot in (0 = new figure each cell)
% tc_width      0.28    fraction of figure width given to the time course panel
% chirp_width   0.56    fraction of figure width given to the chirp panel
%                           (the two panels don't have to sum to 1 - the gap
%                           between them is whatever's left over)
%
% cells are skipped (with a printed note) if either:
%   - the slave cell has no chirp match in datarun1 (see BROWSE_CHIRP), or
%   - the master cell has no time course yet in datarun.stas.time_courses
%     (i.e. GET_STA_SUMMARIES wasn't run for that cell)
%
% See also: LOAD_CHIRP, PLOT_CHIRP_PSTH, PLOT_TIME_COURSE, GET_STA_SUMMARIES,
%           BROWSE_CHIRP, BROWSE_RASTERS, MAP_EI
%
% 2026-07 mcarleton

p = inputParser;
p.addParamValue('figure', 0);
p.addParamValue('tc_width', 0.28);
p.addParamValue('chirp_width', 0.56);
p.parse(varargin{:});
params = p.Results;

if ~isfield(datarun1,'chirp')
    error('browse_chirp_sta: datarun1.chirp not found - run load_chirp(datarun1, chirp_path) first.');
end
if ~isfield(datarun,'stas') || ~isfield(datarun.stas,'time_courses')
    error(['browse_chirp_sta: datarun.stas.time_courses not found - this is not filled in by ' ...
        'load_sta. Run:  datarun = get_sta_summaries(datarun, cell_spec);  first (e.g. with ' ...
        'cell_spec = master_ids_matched).']);
end

master_ids_matched = master_ids_matched(:);
mapped_cell_ids = mapped_cell_ids(:);

if length(master_ids_matched) ~= length(mapped_cell_ids)
    error(['browse_chirp_sta: master_ids_matched (%d) and mapped_cell_ids (%d) are different ' ...
        'lengths - they need to be paired 1:1 (see the header of this function).'], ...
        length(master_ids_matched), length(mapped_cell_ids));
end

% panel layout, normalized figure coordinates [left bottom width height]
% - leaves a margin at the bottom for the status text / Next button
tc_pos    = [0.06, 0.18, params.tc_width,    0.72];
chirp_pos = [0.06 + params.tc_width + 0.06, 0.18, params.chirp_width, 0.72];

n_cells = length(master_ids_matched);
n_skipped = 0;

for k = 1:n_cells
    master_id = master_ids_matched(k);
    slave_id = mapped_cell_ids(k);

    slave_num = get_cell_indices(datarun1, slave_id);
    if ~datarun1.chirp.matched(slave_num)
        fprintf('browse_chirp_sta: skipping master cell %d / slave cell %d (no chirp match)\n', ...
            master_id, slave_id);
        n_skipped = n_skipped + 1;
        continue
    end

    master_num = get_cell_indices(datarun, master_id);
    if isempty(datarun.stas.time_courses{master_num})
        fprintf(['browse_chirp_sta: skipping master cell %d / slave cell %d ' ...
            '(no time course yet - run get_sta_summaries on this cell)\n'], master_id, slave_id);
        n_skipped = n_skipped + 1;
        continue
    end

    % fresh figure with two panels, asymmetric widths
    if isequal(params.figure, 0)
        fig = figure;
    else
        fig = figure(params.figure);
    end
    clf(fig);

    % LEFT (narrow): STA time course. 'clear', false is required here -
    % plot_time_course does a figure-wide clf by default, which would wipe
    % the chirp panel too if left at its default.
    axes('Position', tc_pos); %#ok<LAXES>
    plot_time_course(datarun, master_id, 'figure', -1, 'clear', false);

    % RIGHT (wide): chirp PSTH, phase-shaded (plot_chirp_psth titles itself, incl. QI)
    axes('Position', chirp_pos); %#ok<LAXES>
    plot_chirp_psth(datarun1, slave_id, 'figure', -1);

    % status text + Next button, same pattern as BROWSE_RASTERS/BROWSE_CHIRP
    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02 0.02 0.6 0.05], ...
        'String', sprintf('master %d -> slave %d (%d of %d) - close window to stop', ...
            master_id, slave_id, k, n_cells), ...
        'HorizontalAlignment', 'left', 'FontSize', 10, ...
        'BackgroundColor', get(fig, 'Color'));

    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Next >>', ...
        'Units', 'normalized', 'Position', [0.85 0.02 0.12 0.06], ...
        'FontSize', 11, 'Callback', @(src, evt) uiresume(fig));

    stopRequested = false;
    set(fig, 'CloseRequestFcn', @(src, evt) stopLoop());

    uiwait(fig);

    if stopRequested
        return
    end

    if ishandle(fig)
        close(fig);
    end
end

fprintf('Done browsing %d cell pairs (%d skipped).\n', n_cells, n_skipped);

    function stopLoop()
        stopRequested = true;
        uiresume(fig);
        delete(fig);
    end

end
