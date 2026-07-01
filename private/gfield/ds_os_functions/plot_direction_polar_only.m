function plot_direction_polar_only(tuning_struct, spike_nums, stim_struct, varargin)

p = inputParser;
p.addParameter('fig_num', 1, @isnumeric);
p.addParameter('grating_duration', 8, @isnumeric)
p.addParameter('print', false, @islogical)
p.addParameter('clear_fig', true, @islogical)
p.addParameter('fig_title', []);
p.addParameter('color', 'k')
p.parse(varargin{:});


% ---- BEGIN FUNCTION ---

num_dirs = size(tuning_struct, 1);

directions = stim_struct.params.DIRECTION;
% this is to close the last segment of the polar plot
directions = [directions, directions(1)];
spike_nums = [spike_nums, spike_nums(1)];



% handle case where there are 8 directions
if num_dirs == 8 
    
    
    % plot tuning function
    figure(p.Results.fig_num);
    subplot(3,3,5)
    polarplot(deg2rad(directions), spike_nums);
    
    
    % insert figure title if provided
    if ~isempty(p.Results.fig_title)
        figure(p.Results.fig_num);
        if isstring(p.Results.fig_title)
            title(p.Results.fig_title)
        else
            title(num2str(p.Results.fig_title))
        end
    end
    
    
end

% handle plotting when there are 12 directions
if num_dirs == 12
    
    
    % plot tuning function
    figure(p.Results.fig_num);
    polarplot(deg2rad(directions), spike_nums,'color', p.Results.color);
    drawnow



    % insert figure title if provided
    if ~isempty(p.Results.fig_title)
        figure(p.Results.fig_num);
        if isstring(p.Results.fig_title)
            title(p.Results.fig_title)
        else
            title(num2str(p.Results.fig_title))
        end
    end    
end


end





