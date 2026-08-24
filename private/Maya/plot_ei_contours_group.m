function plot_ei_contours_group(datarun, cell_specification, varargin)
% PLOT_EI_CONTOURS_GROUP   Overlay EI amplitude contours for every cell in a group
%
% usage:  plot_ei_contours_group(datarun, cell_specification, varargin)
%
% arguments:  datarun            - datarun struct with datarun.ei.eis and datarun.ei.position
%                                     already loaded (see LOAD_EI)
%             cell_specification - which cells to plot, see GET_CELL_INDICES for options
%                                     (a type name string, e.g. 'OFFs noise', a cell array of
%                                      type names, a vector of cell ids, or 'all')
%             varargin           - struct or list of optional parameters (see below)
%
% This is a group wrapper around HEX_CONTOUR, which only computes/plots contours for one
% cell at a time. Here each cell in the group is drawn in its own color so you can see how
% a whole cell type tiles the array (e.g. checking mosaic coverage or looking for duplicates).
%
% optional params, their default values, and what they specify:
%
% n_contours        1           number of amplitude contour levels per cell.
%                                   1 = just the outer (soma) outline, which is usually
%                                   clearest when overlaying many cells at once.
% contour_spacing   'linear'    'linear' or 'log' spacing of contour levels (passed to hex_contour)
% colors            []          Nx3 color matrix, one row per cell. if empty, uses HSV
%                                   colormap spaced across the group so every cell is
%                                   distinguishable
% line_width        1.5         width of the plotted contour lines
% plot_electrodes   true        plot electrode positions as small black dots
% array_outline     false       overlay array outline (requires datarun.piece.corners)
% figure            0           figure or axes to plot in, see SET_UP_FIG_OR_AXES
%                                   (0 = new figure, -1 = current axes)
% label             false       label each cell's outer contour with its cell ID
%
% examples:
%
%   plot_ei_contours_group(datarun, 'OFFs noise')
%   plot_ei_contours_group(datarun, {'ON parasol','OFF parasol'}, 'n_contours', 2, 'label', true)
%
% See also: HEX_CONTOUR, PLOT_EI, PLOT_EI_SCROLL, PLOT_RFS_EIS, GET_CELL_INDICES
%
% 2026-07 mcarleton
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addParamValue('n_contours', 1);
p.addParamValue('contour_spacing', 'linear', @(x)any(strcmpi(x,{'linear','log'})));
p.addParamValue('colors', []);
p.addParamValue('line_width', 1.5);
p.addParamValue('plot_electrodes', true);
p.addParamValue('array_outline', false);
p.addParamValue('figure', 0);
p.addParamValue('label', false);
p.parse(varargin{:});
params = p.Results;


% RESOLVE CELLS IN THE GROUP

cell_ids = get_cell_ids(datarun, cell_specification);
n_cells = length(cell_ids);

if n_cells == 0
    warning('plot_ei_contours_group: no cells matched cell_specification.');
    return
end


% SET UP COLORS, ONE PER CELL

if isempty(params.colors)
    colors = hsv(n_cells);
else
    colors = params.colors;
    if size(colors,1) < n_cells
        error('plot_ei_contours_group: need at least one color per cell (%d cells, %d colors given).', ...
            n_cells, size(colors,1));
    end
end


% SET UP AXES

plot_axes = set_up_fig_or_axes(params.figure);
axes(plot_axes);
hold on;

xCoords = datarun.ei.position(:,1);
yCoords = datarun.ei.position(:,2);

if params.plot_electrodes
    plot(xCoords, yCoords, 'ok', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
end

if params.array_outline
    if isfield(datarun,'piece') && isfield(datarun.piece,'corners') && ~isempty(datarun.piece.corners)
        corners = datarun.piece.corners;
        plot(corners([1:end 1],1), corners([1:end 1],2), 'k-');
    else
        warning('plot_ei_contours_group: could not plot array outline, datarun.piece.corners is missing/empty.');
    end
end


% LOOP OVER CELLS, COMPUTE + PLOT CONTOURS

for cc = 1:n_cells

    cell_id = cell_ids(cc);
    ei = get_ei(datarun, cell_id);

    % amplitude used for the contour: max abs value across time, per electrode
    z = max(abs(ei), [], 2);

    % suppress hex_contour's own plotting (fig_or_axes = []); we plot manually below
    % so every cell can get its own color
    contours = hex_contour(xCoords, yCoords, z, params.n_contours, ...
        'contourSpacing', params.contour_spacing, 'fig_or_axes', [], 'plotCoords', false);

    outer_path_for_label = [];

    for lev = 1:length(contours)
        for pp = 1:length(contours(lev).paths)
            path = contours(lev).paths{pp};
            if isempty(path)
                continue
            end
            plot(path(:,1), path(:,2), '-', 'Color', colors(cc,:), 'LineWidth', params.line_width);

            if lev == length(contours) && isempty(outer_path_for_label)
                outer_path_for_label = path;
            end
        end
    end

    if params.label
        if ~isempty(outer_path_for_label)
            label_pos = mean(outer_path_for_label, 1);
        else
            % fall back to the peak electrode if no contour path was found
            [~, peak_elec] = max(z);
            label_pos = [xCoords(peak_elec) yCoords(peak_elec)];
        end
        text(label_pos(1), label_pos(2), num2str(cell_id), 'Color', colors(cc,:), 'FontSize', 8);
    end
end

axis equal
title(sprintf('EI contours: %d cells', n_cells));
hold off
