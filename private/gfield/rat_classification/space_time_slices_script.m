datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

filt_params.radius = 1;

datarun = get_rfs_filtered(datarun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

cell_type_num = get_cell_type_nums(datarun, 'ds');
cell_indices = get_cell_indices(datarun, {cell_type_num});

%%
cell_id = 1098
rgc = find(datarun.cell_types{cell_type_num}.cell_ids == cell_id);

%rgc = 1;
%cell_id = datarun.cell_ids(cell_indices(rgc))



% get the space-time STA
temp_sta = datarun.stas.stas{cell_indices(rgc)};
temp_sta = squeeze(temp_sta);

% apply a spatial filter to STA
if 1
    filt_sta = squeeze(zeros(size(datarun.stas.stas{1})));
    temp_params.radius = 0.75;
    for fm = 1:datarun.stas.depth
        temp_fm = rf_filtered(temp_sta(:,:,fm), temp_params);
        filt_sta(:,:, fm) = temp_fm;
    end
end
filt_sta = reshape(filt_sta, [3200, 30]);

% perform SVD
[U, S, V] = svd(filt_sta);

% calculate and store singular value ratios 
singvals = diag(S);
si = singvals.^2 ./ sum(singvals.^2);
plot(si, 'ko')
axis([0 20 0 1])

% choose number of dimension to recombined to reduce noise
num_dims = 6;
newSTA = U(:,1:num_dims) * S(1:num_dims, 1:num_dims) * V(:, 1:num_dims)';
newSTA = reshape(newSTA, [40 80 30]);

% normalize the new sta
max_sta = max(newSTA(:));
min_sta = min(newSTA(:));
if abs(min_sta) > max(newSTA);
    norm_sta = newSTA ./ abs(min_sta);
else
    norm_sta = newSTA ./ max_sta;
end
norm_sta = (norm_sta +1) ./ 2;

for fm = 1:size(norm_sta,3);
    imagesc(norm_sta(:,:,fm))
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    pause
end



temp_rf = datarun.stas.rfs{cell_indices(1)};
temp_rf = norm_image(temp_rf);
temp_rf = squeeze(temp_rf(:,:,1));
imagesc(temp_rf)
colormap(brewermap([],'RdBu'))
caxis([0,1]) 

%get pixel at center of mass
temp_rf_com = round(datarun.stas.rf_coms{cell_indices(1)});

temp_sta = datarun.stas.stas{cell_indices(1)};
slice_1 = squeeze(temp_sta(temp_rf_com(2)-10:temp_rf_com(2)+10,temp_rf_com(1),:,:));
slice_2 = squeeze(temp_sta(temp_rf_com(2), temp_rf_com(1)-10:temp_rf_com(1)+10,:,:));

st_plot = (slice_1 + slice_2)/2;

st_plot = norm_image(st_plot);
st_plot = squeeze(st_plot(:,:,1));
imagesc(st_plot)
colormap(brewermap([],'RdBu'))
caxis([0,1]) 
