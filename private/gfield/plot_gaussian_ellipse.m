function h = plot_gaussian_ellipse(fitParams, varargin)
% PLOT_GAUSSIAN_ELLIPSE Plot ellipse(s) representing 2D Gaussian fit contours
%
% Inputs:
%   fitParams - Structure from fit_2d_gaussian_rf containing:
%               .x0, .y0 - center coordinates
%               .sigma_x, .sigma_y - standard deviations
%               .theta - rotation angle
%
% Optional Name-Value Pairs:
%   'SigmaLevel' - Sigma level(s) for ellipse (default: 1)
%                  Can be a scalar or vector for multiple ellipses
%   'Color' - Line color (default: 'r')
%   'LineWidth' - Line width (default: 2)
%   'LineStyle' - Line style (default: '-')
%   'ShowCenter' - Show center marker (default: true)
%   'CenterMarker' - Marker style for center (default: '+')
%   'CenterSize' - Marker size for center (default: 15)
%   'NumPoints' - Number of points for ellipse (default: 100)
%   'Labels' - Show sigma level labels (default: false)
%   'LabelOffset' - Offset for labels in pixels (default: 5)
%
% Output:
%   h - Structure containing handles to plotted objects:
%       .ellipses - handles to ellipse lines
%       .center - handle to center marker
%       .labels - handles to text labels (if Labels is true)
%
% Examples:
%   % Plot 1-sigma ellipse
%   plot_gaussian_ellipse(fitParams);
%
%   % Plot multiple sigma levels
%   plot_gaussian_ellipse(fitParams, 'SigmaLevel', [1, 2, 3]);
%
%   % Customize appearance
%   plot_gaussian_ellipse(fitParams, 'SigmaLevel', 2, ...
%                         'Color', 'b', 'LineWidth', 3, ...
%                         'ShowCenter', true, 'Labels', true);

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'fitParams');
    addParameter(p, 'SigmaLevel', 1, @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'Color', 'k');
    addParameter(p, 'LineWidth', 1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'LineStyle', '-');
    addParameter(p, 'ShowCenter', false, @islogical);
    addParameter(p, 'CenterMarker', 'o');
    addParameter(p, 'CenterSize', 15, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'NumPoints', 100, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Labels', false, @islogical);
    addParameter(p, 'LabelOffset', 5, @(x) isnumeric(x));
    
    parse(p, fitParams, varargin{:});
    opts = p.Results;
    
    % Extract parameters
    x0 = fitParams.x0;
    y0 = fitParams.y0;
    sigma_x = fitParams.sigma_x;
    sigma_y = fitParams.sigma_y;
    theta = fitParams.theta;
    
    % Initialize output structure
    h.ellipses = [];
    h.center = [];
    h.labels = [];
    
    % Hold current plot
    wasHeld = ishold;
    hold on;
    
    % Generate parametric ellipse
    t = linspace(0, 2*pi, opts.NumPoints);
    
    % Plot ellipse for each sigma level
    nLevels = length(opts.SigmaLevel);
    colors = get_color_gradient(opts.Color, nLevels);
    
    for i = 1:nLevels
        sigma_level = opts.SigmaLevel(i);
        
        % Parametric ellipse at this sigma level
        ellipse_x = sigma_level * sigma_x * cos(t);
        ellipse_y = sigma_level * sigma_y * sin(t);
        
        % Rotate ellipse
        rot_x = ellipse_x * cos(theta) - ellipse_y * sin(theta);
        rot_y = ellipse_x * sin(theta) + ellipse_y * cos(theta);
        
        % Translate to center
        x_plot = x0 + rot_x;
        y_plot = y0 + rot_y;
        
        % Determine color for this ellipse
        if nLevels > 1
            current_color = colors(i, :);
        else
            current_color = opts.Color;
        end
        
        % Plot ellipse
        h.ellipses(i) = plot(x_plot, y_plot, ...
            'Color', current_color, ...
            'LineWidth', opts.LineWidth, ...
            'LineStyle', opts.LineStyle);
        
        % Add label if requested
        if opts.Labels
            % Find point on ellipse for label (at angle pi/4)
            label_angle = pi/4;
            label_x_ellipse = sigma_level * sigma_x * cos(label_angle);
            label_y_ellipse = sigma_level * sigma_y * sin(label_angle);
            
            % Rotate
            label_x_rot = label_x_ellipse * cos(theta) - label_y_ellipse * sin(theta);
            label_y_rot = label_x_ellipse * sin(theta) + label_y_ellipse * cos(theta);
            
            % Translate
            label_x = x0 + label_x_rot;
            label_y = y0 + label_y_rot;
            
            % Create label text
            if sigma_level == 1
                label_text = '1\sigma';
            else
                label_text = sprintf('%.1f\\sigma', sigma_level);
            end
            
            h.labels(i) = text(label_x, label_y, label_text, ...
                'Color', current_color, ...
                'FontWeight', 'bold', ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end
    end
    
    % Plot center marker
    if opts.ShowCenter
        h.center = plot(x0, y0, opts.CenterMarker, ...
            'Color', opts.Color, ...
            'MarkerSize', opts.CenterSize, ...
            'LineWidth', opts.LineWidth);
    end
    
    % Restore hold state
    if ~wasHeld
        hold off;
    end
end

function colors = get_color_gradient(base_color, n)
% Generate color gradient for multiple ellipses
% Lighter colors for larger sigma levels
    
    if ischar(base_color) || isstring(base_color)
        % Convert color name to RGB
        switch lower(base_color)
            case 'r'
                base_rgb = [1, 0, 0];
            case 'g'
                base_rgb = [0, 1, 0];
            case 'b'
                base_rgb = [0, 0, 1];
            case 'c'
                base_rgb = [0, 1, 1];
            case 'm'
                base_rgb = [1, 0, 1];
            case 'y'
                base_rgb = [1, 1, 0];
            case 'k'
                base_rgb = [0, 0, 0];
            case 'w'
                base_rgb = [1, 1, 1];
            otherwise
                base_rgb = [1, 0, 0]; % default to red
        end
    else
        base_rgb = base_color;
    end
    
    if n == 1
        colors = base_rgb;
    else
        % Create gradient from base color to lighter version
        colors = zeros(n, 3);
        for i = 1:n
            % Darker for inner ellipses, lighter for outer
            alpha = 1 - (i-1)/(n-1) * 0.6; % Range from 1.0 to 0.4
            colors(i, :) = base_rgb * alpha + [1, 1, 1] * (1 - alpha);
        end
    end
end
