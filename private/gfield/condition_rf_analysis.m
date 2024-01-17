datapath = '/Users/gfield/Analysis/rat/2012-10-31-0/data000-map/data000-map';
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/Chichilnisky-lab/2017-11-29-0/data000/data000';
datapath = '/Users/gfield/Analysis/rat/2012-10-15-0/data000-map/data000-map';
datapath = '/Users/gfield/Analysis/rat/2012-10-10-1/data001-3600-7200s/data001-map/data001-map';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = set_polarities(datarun);

marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);


%filt_rad = 0.75;
%datarun = get_rfs_filtered(datarun,'all','verbose',1,'filt_params',struct('filt_type','gauss','radius',filt_rad));


%%

on_bt_indices = get_cell_indices(datarun, 'OFF type2');
off_bt_indices = get_cell_indices(datarun, 'ON type2');

on_bt_COMs = zeros(length(on_bt_indices), 2);
for rgc = 1:length(on_bt_indices)
    on_bt_COMs(rgc, :) = datarun.stas.fits{on_bt_indices(rgc)}.mean;
end

off_bt_COMs = zeros(length(off_bt_indices), 2);
for rgc = 1:length(off_bt_indices)
    off_bt_COMs(rgc, :) = datarun.stas.fits{off_bt_indices(rgc)}.mean;
end


num_neighbors = length(off_bt_indices);
%num_neighbors = 10;
up_scale_factor = 10;
polarity = -1;
[rf_height, rf_width] = size(datarun.stas.rfs{1});
rf_height = rf_height * up_scale_factor;
rf_width = rf_width * up_scale_factor;
window_size = 6 * up_scale_factor;

summed_rfs = zeros(2*window_size+1, 2*window_size+1);

for rgc = 1:length(on_bt_indices)
    
    d = ipdm(on_bt_COMs(rgc, :), off_bt_COMs);
    [sort_d, tmp_inds] = sort(d, 'ascend');
    
    on_bt_id = datarun.cell_ids(on_bt_indices(rgc));
    off_bt_ids = datarun.cell_ids(off_bt_indices(tmp_inds(1:num_neighbors)));
    
    figure(1); clf;
    plot_rf_summaries(datarun, on_bt_id, 'plot_fits', true, 'fit_color', 'b', 'foa', 1, 'clear', true)
    plot_rf_summaries(datarun, off_bt_ids, 'plot_fits', false, 'fit_color', 'r', 'foa', 1, 'clear', false)
    hold on
    plot(on_bt_COMs(rgc,1), on_bt_COMs(rgc,2), 'bo')
    [on_bt_COMs(rgc,1), on_bt_COMs(rgc,2)]
    
    user_input = input('do you want to keep this group for analysis?');
    
    if user_input == 1
    
        % compute window
        tmp_origin = round(on_bt_COMs(rgc,:) * up_scale_factor);
        
        bg_x = tmp_origin(1) - window_size;
        ed_x = tmp_origin(1) + window_size;
        bg_y = tmp_origin(2) - window_size;
        ed_y = tmp_origin(2) + window_size;
        
        if bg_x < 1; bg_x = 1; end
        if ed_x > rf_width; ed_x = rf_width; end
        if bg_y < 1; bg_y = 1; end
        if ed_y > rf_height; ed_y = rf_height; end
     
        for cc = 1:num_neighbors
            
            tmp_rf = datarun.stas.rfs{off_bt_indices(tmp_inds(cc))};
            % tmp_rf = datarun.stas.summaries_filt{off_bt_indices(tmp_inds(cc))};
            if length(size(tmp_rf)) == 3
                tmp_rf = squeeze(tmp_rf(:,:,2));
            end
            norm_rf = reshape(tmp_rf, 1, []);
            norm_rf = norm_rf ./ (polarity * ext(norm_rf));
            norm_rf = norm_rf ./ std(norm_rf);
            norm_rf = reshape(norm_rf, size(tmp_rf));
            norm_rf = matrix_scaled_up(norm_rf, up_scale_factor);
            norm_rf = norm_rf(bg_y:ed_y,bg_x:ed_x);
            
%             figure(1); clf
%             imagesc(norm_rf)
%             pause
            
            summed_rfs = summed_rfs + norm_rf;
            figure(2); clf;
            imagesc(summed_rfs)
            pause(0.1)
        end
    else
        continue
        
    end
end
    
    
    
%% get radial average
distance_list = zeros(size(summed_rfs));

cntr_pt = 63;

for x_ind = 1:size(distance_list,1)
    for y_ind = 1:size(distance_list,2)        
        distance_list(x_ind, y_ind) = sqrt((x_ind-cntr_pt).^2 + (y_ind-cntr_pt).^2);
    end
end

[sorted_dists, sort_inds] = sort(distance_list(:), 'ascend');
figure(4)
plot(sorted_dists,summed_rfs(sort_inds), '.')

ordered_energies = summed_rfs(sort_inds);

partition_factor = 20;
dist_ranges = 0:(cntr_pt./partition_factor):cntr_pt;
rad_energy_mean = zeros(length(dist_ranges)-1,1);
rad_energy_ste = zeros(length(dist_ranges)-1,1);
for cnt = 1:(length(dist_ranges)-1)
    temp_inds = find(sorted_dists > dist_ranges(cnt) & sorted_dists < dist_ranges(cnt+1));
    rad_energy_mean(cnt) = mean(ordered_energies(temp_inds));
    rad_energy_ste(cnt) = std(ordered_energies(temp_inds)) ./ sqrt(length(temp_inds));
end

x_units = 1/partition_factor:1/partition_factor:1;
figure(5); clf;
errorbar(x_units, rad_energy_mean, rad_energy_ste, 'ko-')
hold on
plot([0 2], [0 0], 'k')
ylabel('delta zscore')
xlabel('normalized distance')
hold off

    
    
    
    