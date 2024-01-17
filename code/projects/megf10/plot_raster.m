function plot_raster(raster, start_time, end_time, varargin)
% plot_raster(raster, start_time, end_time)
% raster: nx1 cells, one cell for a trial
%
% xyao
% 2013-12-16

p = inputParser;
p.addParameter('raster_color', 'k', @isstring)
p.parse(varargin{:});

    
for j = 1:length(raster)
    SpikeTime = raster{j};
    SpikeTime = SpikeTime';
    X = [SpikeTime; SpikeTime];
    Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
    line(X, Y, 'color', p.Results.raster_color);
    axis([start_time, end_time,0,length(raster)])
    hold on
end
