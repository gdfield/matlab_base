function browse_contrast_response(crf, varargin)
%BROWSE_CONTRAST_RESPONSE  Step through cells one at a time, plotting each
%contrast-response function with PLOT_CONTRAST_RESPONSE and advancing to
%the next cell when you click "Next".
%
%   browse_contrast_response(crf)
%
% INPUTS
%   crf       - struct array from GET_CONTRAST_RESPONSE (one entry per cell)
%   varargin  - any name-value pairs accepted by PLOT_CONTRAST_RESPONSE,
%               passed straight through for every cell
%
% BEHAVIOR: same click-through pattern as BROWSE_RASTERS / BROWSE_CHIRP -
% plots the first cell, shows a "Next >>" button, closing the window stops
% the loop early.
%
% See also: GET_CONTRAST_RESPONSE, PLOT_CONTRAST_RESPONSE, BROWSE_RASTERS, BROWSE_CHIRP
%
% 2026-07 mcarleton

n_cells = length(crf);

for i = 1:n_cells
    plot_contrast_response(crf, 'cell_index', i, varargin{:});
    fig = gcf;

    % status text so you know where you are in the list
    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02 0.02 0.5 0.05], ...
        'String', sprintf('Cell %d (%d of %d) - close window to stop', ...
            crf(i).cell_id, i, n_cells), ...
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

fprintf('Done browsing %d contrast-response functions.\n', n_cells);

    function stopLoop()
        stopRequested = true;
        uiresume(fig);
        delete(fig);
    end

end
