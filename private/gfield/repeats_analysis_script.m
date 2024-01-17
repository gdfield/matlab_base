base_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-02-15-0/';
ndf{1} = [base_path, 'data002-map/data002-map'];
ndf{2} = [base_path, 'data004-map/data004-map'];
ndf{3} = [base_path, 'data006-map/data006-map'];
ndf{4} = [base_path, 'data008-map/data008-map'];
ndf{5} = [base_path, 'data009-map/data009-map'];

num_light_levels = length(ndf);

for ll = 1:num_light_levels
    temp_datarun = load_data(ndf{ll});
    temp_datarun = load_neurons(temp_datarun);
    temp_datarun = load_ei(temp_datarun, 'all');
    datarun{ll} = temp_datarun;
end

master_path = [base_path, 'data010/data010'];
master_datarun = load_data(master_path);
master_datarun = load_neurons(master_datarun);
master_datarun = load_ei(master_datarun, 'all');
master_datarun = load_params(master_datarun);


cell_type = 'off brisk transient';
temp_cell_indices = get_cell_indices(master_datarun, cell_type)

for ll = 1:num_light_levels    
    temp_list = map_ei(master_datarun, datarun{ll}, 'master_cell_type', cell_type);  
    cell_list{ll} = temp_list(temp_cell_indices);
end



cd ~/Desktop/2019-02-15/

for rgc = 1:length(temp_cell_indices)

    figure(1); clf;
    
    for ll = 1:num_light_levels
        %trig_indices = find(diff(datarun{ll}.triggers) > 1.69);    
        trig_indices = 1:6:300;
        
        if ~isempty(cell_list{ll}{rgc})
            temp_index = get_cell_indices(datarun{ll}, cell_list{ll}{rgc});
            temp_spikes = datarun{ll}.spikes{temp_index};   
            subplot(5, 1, ll)
            raster_struct = get_raster(temp_spikes, datarun{ll}.triggers(trig_indices(1:50)), 'stop', 9, 'foa', 1, 'labels', 0);
            axis([0 9 0 50])
            yticks([0 25 50])
            ylabel('trials')
        end
    end
    subplot(5,1,1)
    temp_title = ['Cell-' num2str(master_datarun.cell_ids(temp_cell_indices(rgc)))];
    title(temp_title)
    subplot(5,1,5)
    xlabel('seconds')
    drawnow
    %reply = input('export figure? [y/n]', 's');
    %if strcmp(reply, 'y') 
        temp_filename = [temp_title, '.pdf']
        exportgraphics(gcf, temp_filename, 'ContentType', 'vector')
    %end
end

    
    


