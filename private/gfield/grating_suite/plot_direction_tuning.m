function plot_direction_tuning(tuning_struct, spike_nums, stim_struct, varargin)

p = inputParser;
p.addParameter('fig_num', 1, @isnumeric);
p.addParameter('grating_duration', 8, @isnumeric)
p.addParameter('print', false, @islogical)
p.addParameter('clear_fig', true, @islogical)
p.addParameter('save_path', '~/Desktop/', @ischar)
p.addParameter('save_name', 'tuning_plot.pdf', @ischar)
p.addParameter('fig_title', []);
p.addParameter('print_for_fig', false, @islogical)
p.parse(varargin{:});


% ---- BEGIN FUNCTION ---

num_dirs = size(tuning_struct, 1);
num_repeats = size(tuning_struct, 2);
g_duration = p.Results.grating_duration;

directions = stim_struct.params.DIRECTION;
% this is to close the last segment of the polar plot
directions = [directions, directions(1)];
spike_nums = [spike_nums, spike_nums(1)];


figure(p.Results.fig_num)
if p.Results.clear_fig
    clf;
end

% handle case where there are 8 directions
if num_dirs == 8 
    
    % this alignes the raster locations with the polar plot
    subplot_ind = [6, 3, 2, 1, 4, 7, 8, 9];
    sep_plot_nums = 0:45:315;
    
    % plot tuning function
    figure(p.Results.fig_num);
    subplot(3,3,5)
    polarplot(deg2rad(directions), spike_nums);
    
    for g_dir = 1:num_dirs
        subplot(3,3,subplot_ind(g_dir))
        tmp_spike_times = tuning_struct(g_dir, :);
        plot_raster(tmp_spike_times, 0, g_duration)
        axis square
        xlim([0 g_duration])
    end
    
    % insert figure title if provided
    if ~isempty(p.Results.fig_title)
        figure(p.Results.fig_num);
        subplot(3,3,1)
        if isstring(p.Results.fig_title)
            title(p.Results.fig_title)
        else
            title(num2str(p.Results.fig_title))
        end
    end
    
    
end

% handle plotting when there are 12 directions
if num_dirs == 12
    
    subplot_vector = [13 15 10 4 3 2 6 11 16 22 23 24 20]; 
    sep_plot_nums = 0:30:330;
    
    % plot tuning function
    figure(p.Results.fig_num);
    subplot(3,3,5)
    polarplot(deg2rad(directions), spike_nums);
    drawnow

    if p.Results.print_for_fig
        figure(50)
        polarplot(deg2rad(directions), spike_nums);
        print(50, '~/Desktop/polar.pdf', '-dpdf')
    end

    for g_dir = 1:num_dirs
        figure(p.Results.fig_num);
        subplot(5,5,subplot_vector(g_dir+1))
        tmp_spike_times = tuning_struct(g_dir, :);
        plot_raster(tmp_spike_times, 0, g_duration)
        axis square
        xlim([0 g_duration])
        if g_dir == 1
            set(gca,'xtick',[0 g_duration])
            set(gca,'ytick',[0 size(tuning_struct, 2)])
        else
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'xticklabel',[])
            set(gca,'yticklabel',[])
        end

        
        if p.Results.print_for_fig
            figure(10+g_dir)
            plot_raster(tmp_spike_times, 0, g_duration)
            figure_title = ['direction ', num2str(sep_plot_nums(g_dir))];
            title(figure_title)
            save_name = ['~/Desktop/direction', num2str(sep_plot_nums(g_dir)),'.pdf'];
            print(10+g_dir, save_name, '-dpdf')
        end
        
    end
    % insert figure title if provided
    if ~isempty(p.Results.fig_title)
        figure(p.Results.fig_num);
        subplot(3,3,5)
        if isstring(p.Results.fig_title)
            title(p.Results.fig_title)
        else
            title(num2str(p.Results.fig_title))
        end
    end    
end


if p.Results.print 
    % construct save location and filename
    save_final = [p.Results.save_path, p.Results.save_name];
    %print(p.Results.fig_num, save_final, '-dpdf', '-bestfit')
    
end





