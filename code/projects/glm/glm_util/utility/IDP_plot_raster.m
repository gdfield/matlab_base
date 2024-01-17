function  IDP_plot_raster(prepped_data, cell_to_plot)

% prepped data must be the output of interleaved_data_prep
% cell_to_plot is the index of the cell in the prepped data structure

figure;
set(gcf, 'Position', 100*[1 1 10 3])
hold on
for i = 1:size(prepped_data.testspikes, 1)
    plot(prepped_data.testspikes{i, cell_to_plot}, i, 'k.')
end

end