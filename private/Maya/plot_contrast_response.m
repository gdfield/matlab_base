function plot_contrast_response(crf, varargin)
% PLOT_CONTRAST_RESPONSE   Plot mean firing rate vs. contrast, with SEM error bars
%
% usage: plot_contrast_response(crf, varargin)
%
% arguments: crf      - struct array from GET_CONTRAST_RESPONSE (one entry per cell)
%            varargin - struct or list of optional parameters (see below)
%
% optional params, their default values, and what they specify:
%
% figure       0     figure or axes to plot in, see SET_UP_FIG_OR_AXES
%                        (0 = new figure, -1 = current axes)
% cell_index   []    which entry/entries of crf to plot, by position (default:
%                        all of them, overlaid with a legend)
%
% See also: GET_CONTRAST_RESPONSE, BROWSE_CONTRAST_RESPONSE
%
% 2026-07 mcarleton
%

p = inputParser;
p.addParamValue('figure', 0);
p.addParamValue('cell_index', []);
p.parse(varargin{:});
params = p.Results;

if isempty(params.cell_index)
    idx = 1:length(crf);
else
    idx = params.cell_index;
end

plot_axes = set_up_fig_or_axes(params.figure);
axes(plot_axes);
hold on;

colors = hsv(max(length(idx),1));
legend_labels = cell(1, length(idx));

for ii = 1:length(idx)
    c = crf(idx(ii));
    errorbar(c.contrast, c.mean_rate, c.sem_rate, '-o', 'Color', colors(ii,:), ...
        'MarkerFaceColor', colors(ii,:), 'LineWidth', 1.5);
    legend_labels{ii} = sprintf('cell %d', c.cell_id);
end

xlabel('contrast');
ylabel('firing rate (Hz)');
if length(idx) == 1
    title(sprintf('contrast response - cell %d', crf(idx(1)).cell_id));
else
    title('contrast response function');
    legend(legend_labels, 'Location', 'best');
end
hold off;
