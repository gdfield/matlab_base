function plot_raster(raster, start_time, end_time, varargin)
% plot_raster(raster, start_time, end_time)
% raster: nx1 cells, one cell for a trial
%
% xyao
% 2013-12-16

p = inputParser;
p.addParameter('raster_color', 'k')
p.addParameter('fast', false', @islogical)
p.parse(varargin{:});
   

if p.Results.fast
    num_trials = length(raster);
    for tr = 1:num_trials
        if ~isempty(raster{tr})
            plot(raster{tr},tr, 'o','MarkerEdgeColor',p.Results.raster_color)
        end
        hold on
    end
    axis([start_time, end_time, 0 num_trials])
else

    for j = 1:length(raster)
        SpikeTime = raster{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', p.Results.raster_color);
        axis([start_time, end_time,0,length(raster)])
        hold on
    end

end
