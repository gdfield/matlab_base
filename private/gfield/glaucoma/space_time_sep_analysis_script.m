%datapath = '/Volumes/FieldLab/Chichilnisky/Analysis/controls/20130401A/chunk1/kilosort25/kilosort25';


datapaths = get_nerve_crush_data('cohort', 'control');

num_dsets = length(datapaths.paths);
for dset = 1:num_dsets

    % load data
    datarun = load_data(datapaths.paths{dset});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    datarun = load_sta(datarun, 'load_sta', 'all');
    datarun = load_ei(datarun, 'all', 'array_type', 512);
    
    % setup some stuff
    num_rgcs = length(datarun.cell_ids);
    marks_params.thresh = 3.5;
    datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);
    
    % set filter radius
    filt_rad = 0.75; % units of stixels
    verbose = 1;
    num_dims_to_plot = 4;
    num_time_steps = 24;
    starting_frame = 30 - num_time_steps;
    microns_per_pixel = 5.5;
    pixels_per_stixel = datarun.stimulus.stixel_height;
    
    
    cell_type = 'all';
    temp_indices = get_cell_indices(datarun, cell_type);
    
    num_rgcs = length(temp_indices);
    cnt = 1;
    % loop through the STAs for each cell and compute the SVD
    for rgc = 1:num_rgcs
        temp_sta = datarun.stas.stas{temp_indices(rgc)};
    
        [x_dim, y_dim, color_dim, num_frames]= size(temp_sta);
    
        % check the color dimension and squeeze accordinging
        if color_dim == 1
            temp_sta = squeeze(temp_sta);
    
        elseif color_dim == 3
            temp_sta = squeeze(temp_sta(:,:,2,:));
    
        end
    
        filt_sta = zeros(x_dim, y_dim, num_time_steps);
        for fm = 1:num_time_steps
            temp_frame = rf_filtered(temp_sta(:,:,fm+starting_frame), 'radius', filt_rad);
            filt_sta(:,:,fm) = temp_frame;
        end
        
        % Perform SVD on the filtered STAs
        [U, S, V] = svd(reshape(filt_sta, [], num_time_steps));
        
        % Store the SVD results in the datarun structure
        datarun.stas.svd{rgc} = struct('U', U(:,1:num_dims_to_plot), 'S', diag(S), 'V', V(:,1:num_dims_to_plot));
    
        % Compute the explained variance for each RGC
        explainedVariance = diag(S).^2 / sum(diag(S).^2);
        datarun.stas.explained_variance{rgc} = explainedVariance;
    
        if verbose
    
            % view the first several dimensions of the SVD
            figure(1)
            for pl = 1:num_dims_to_plot
                subplot(2,4,2*pl-1)
                imagesc(reshape(U(:,pl), [x_dim, y_dim]));
                title(num2str(pl))
                axis off
                axis square
                subplot(2, 4, 2*pl)
                plot(V(:,pl))
                axis off
                axis square
            end
            
            % view the singular vals in terms of var explained
            figure(2)
            plot(explainedVariance, 'o')
            temp_title = ['cell ID is ', num2str(datarun.cell_ids(temp_indices(rgc)))];
            title(temp_title)
            xlabel('dimensions')
            ylabel('variance explained')
        
            pause
        end
    
        if explainedVariance(1) > 0.67
            separable_cell_list(cnt) = temp_indices(rgc);
            cnt = cnt +1;
        end
    
    end
    
    sep_cell_ids = datarun.cell_ids(separable_cell_list);
    num_sep_rgcs = length(separable_cell_list);
    
    
    rf_radii = get_rf_fit_radius(datarun, sep_cell_ids, 'fits_to_use', 'vision', 'units', 'microns');
    rf_areas = get_rf_areas_from_marks(datarun, sep_cell_ids);
    
    
    ttz = zeros(1,num_sep_rgcs);
    for rgc = 1:num_sep_rgcs
    
        temp_tc = datarun.stas.time_courses{separable_cell_list(rgc)};
    
        [tc_fit, final_params] = fit_time_course(temp_tc);
            time_range = 4:0.01:10;
            zoomed_fit = fit_time_course_function(final_params(1:6), time_range);
            [~, temp_ttz] = min(abs(zoomed_fit));
            %ttz(rgc) = ((temp_ttz * 0.01)+4) * (p.Results.ms_per_frame);
            ttz(rgc) = ((temp_ttz * 0.01)+4) * 16.7;
    end
    
    rf_summaries(dset).rf_radii = rf_radii;
    rf_summaries(dset).rf_areas = rf_areas;
    rf_summaries(dset).ttz = ttz;

end

for dset = 1:num_dsets
    % histogram of times-to-zero
    % cumulative
    figure(1)
    [temp_counts, bin_edges] = histcounts(rf_summaries(dset).ttz, 50);
    ttz_cumulative = cumsum(temp_counts) ./ sum(temp_counts);
    plot(bin_edges(2:51), ttz_cumulative)
    hold on
    title('time to zero')
    
    % histogram of RF radii
    figure(2)
    [temp_counts, bin_edges] = histcounts(rf_summaries(dset).rf_radii, 50);
    ttz_cumulative = cumsum(temp_counts) ./ sum(temp_counts);
    plot(bin_edges(2:51), ttz_cumulative)
    hold on
    title('rf sizes')
    
end
