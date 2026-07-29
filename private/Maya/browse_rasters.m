function browse_rasters(datarun, cell_ids, varargin)
%BROWSE_RASTERS  Step through cells one at a time, plotting each with
%TESTRASTER and advancing to the next cell when you click "Next".
%
%   browse_rasters(datarun)
%   browse_rasters(datarun, cell_ids)
%   browse_rasters(datarun, cell_ids, 'direction', 90)
%   browse_rasters(datarun, cell_ids, 'direction', 90, 'temporal_period', 256)
%
% INPUTS
%   datarun    - the datarun struct (after load_data/load_neurons/
%                load_params/load_stim)
%   cell_ids   - vector of cell IDs to step through, in order (default:
%                datarun.cell_ids, i.e. every cell)
%   varargin   - any name-value pairs accepted by TESTRASTER (e.g.
%                'direction', 'temporal_period', 'bin_width'), passed
%                straight through for every cell
%
% BEHAVIOR
%   Plots the first cell using TESTRASTER, then shows a "Next >>" button
%   on the figure. Clicking it closes the current figure and plots the
%   next cell in the list. Closing the figure window directly (the X
%   button) stops the loop early.

if nargin < 2 || isempty(cell_ids)
    cell_ids = datarun.cell_ids;
end

for i = 1:length(cell_ids)
    cell_id = cell_ids(i);

    testraster(datarun, cell_id, varargin{:});
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

fprintf('Done browsing all %d cells.\n', length(cell_ids));

    function stopLoop()
        stopRequested = true;
        uiresume(fig);
        delete(fig);
    end

end
