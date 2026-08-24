function browse_chirp(datarun, cell_ids, varargin)
%BROWSE_CHIRP  Step through cells one at a time, plotting each with
%PLOT_CHIRP_PSTH and advancing to the next cell when you click "Next".
%
%   browse_chirp(datarun)
%   browse_chirp(datarun, cell_ids)
%   browse_chirp(datarun1, mapped_cell_ids)
%
% INPUTS
%   datarun    - datarun struct with datarun.chirp already populated (see
%                LOAD_CHIRP). Note: if you're passing mapped_cell_ids from
%                MAP_EI, make sure this is whichever datarun the chirp file
%                was actually loaded onto (usually the slave/mapped datarun,
%                e.g. datarun1, not the WN master) - mapped_cell_ids are IDs
%                in that datarun's ID space.
%   cell_ids   - vector of cell IDs to step through, in order (default:
%                datarun.cell_ids, i.e. every cell)
%   varargin   - any name-value pairs accepted by PLOT_CHIRP_PSTH, passed
%                straight through for every cell
%
% cells with no chirp match (datarun.chirp.matched == false, e.g. because
% they weren't in the chirp file, or mapped_cell_ids has an ID this datarun's
% chirp data doesn't cover) are skipped with a printed note, rather than
% erroring out the whole browse.
%
% BEHAVIOR: same click-through pattern as BROWSE_RASTERS - plots the first
% matched cell, shows a "Next >>" button, closing the window stops early.
%
% See also: LOAD_CHIRP, PLOT_CHIRP_PSTH, PLOT_CHIRP_RASTER, BROWSE_RASTERS
%
% 2026-07 mcarleton

if ~isfield(datarun,'chirp')
    error('browse_chirp: datarun.chirp not found - run load_chirp(datarun, chirp_path) first.');
end

if nargin < 2 || isempty(cell_ids)
    cell_ids = datarun.cell_ids;
end

n_skipped = 0;

for i = 1:length(cell_ids)
    cell_id = cell_ids(i);

    cell_num = get_cell_indices(datarun, cell_id);
    if ~datarun.chirp.matched(cell_num)
        fprintf('browse_chirp: skipping cell %d (no chirp match in this datarun)\n', cell_id);
        n_skipped = n_skipped + 1;
        continue
    end

    plot_chirp_psth(datarun, cell_id, varargin{:});
    fig = gcf;

    % status text so you know where you are in the list
    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02 0.02 0.5 0.05], ...
        'String', sprintf('Cell %d (%d of %d) - close window to stop', ...
            cell_id, i, length(cell_ids)), ...
        'HorizontalAlignment', 'left', 'FontSize', 10, ...
        'BackgroundColor', get(fig, 'Color'));

    % "Next" button - resumes execution when clicked
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Next >>', ...
        'Units', 'normalized', 'Position', [0.85 0.02 0.12 0.06], ...
        'FontSize', 11, 'Callback', @(src, evt) uiresume(fig));

    % if the user closes the window instead of clicking Next, stop the loop
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

fprintf('Done browsing chirp PSTHs (%d of %d skipped, no chirp match).\n', n_skipped, length(cell_ids));

    function stopLoop()
        stopRequested = true;
        uiresume(fig);
        delete(fig);
    end

end
