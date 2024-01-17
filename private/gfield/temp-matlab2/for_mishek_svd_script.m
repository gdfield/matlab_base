datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2012-10-15-0/data000/data000';

datarun = load_data(datapath);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_neurons(datarun);

cell_types = {'OFF brisk transient'}
num_types = length(cell_types);
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, cell_types, 'marks_params', marks_params, 'sta', 'stas');

cell_indices = get_cell_indices(datarun, cell_types);

filt_params.radius = 0.75;
frames_to_use = 20;
fm_window = 14;
[height, width, colors, nframes] = size(datarun.stas.stas{cell_indices(1)});
verbose = true;

SVDVarExpl = [];
n_dims = 3;
for cts = 1:num_types
    cell_indices = get_cell_indices(datarun, cell_types{cts});
    temp_num_cells = length(cell_indices);
    
    tmp_var_exp = zeros(temp_num_cells,1);
    for rgc = 1:temp_num_cells
        
        %get sta for cell and squeeze
        tmp_sta = datarun.stas.stas{cell_indices(rgc)};
        tmp_sta = squeeze(tmp_sta(:,:,:,nframes - frames_to_use + 1:nframes));
        
        %filter frames
        for fm = 1:size(tmp_sta,3)
            tmp_sta(:,:,fm) = rf_filtered(squeeze(tmp_sta(:,:,fm)), filt_params);
        end
        
        %cut out a portion of the STA
        rf_com = round(datarun.vision.sta_fits{cell_indices(rgc)}.mean);
        h_begin = (height-rf_com(2)) - fm_window;
        h_end = (height-rf_com(2)) + fm_window;
        w_begin = rf_com(1) - fm_window;
        w_end = rf_com(1) + fm_window;
        if h_begin < 1
            h_begin = 1;
            h_end = fm_window*2 + 1;
        end
        if w_begin < 1
            w_begin = 1;
            w_end = fm_window*2 + 1;
        end
        if h_end > height
            h_end = height;
            h_begin = height - 2*fm_window;
        end
        if w_end > width
            w_end = width;
            w_begin = width - 2*fm_window;
        end
        
        tmp_sta = tmp_sta(h_begin:h_end, w_begin:w_end,:);
        new_sizes = size(tmp_sta);

        %reshape and perfrom SVD
        tmp_sta = reshape(tmp_sta, [new_sizes(1)*new_sizes(2), frames_to_use]);        
        [U, S, V] = svd(tmp_sta);

        %calculate and store singular value ratios as explained variance 
        singvals = diag(S);
        si = singvals.^2 ./ sum(singvals.^2);
        
        if verbose
            figure(10)
            semilogy(si, 'ko')
            axis([0 20 0 1])

            for nd = 1:n_dims
                %normalize the new sta
                max_sta = max(U(:,nd));
                min_sta = min(U(:,nd));
                if abs(min_sta) > max_sta
                    norm_rf = U(:,nd) ./ abs(min_sta);
                else
                    norm_rf = U(:,nd) ./ max_sta;
                end
                norm_rf = (norm_rf +1) ./ 2;


                figure(11)
                subplot(3,1,nd)
                imagesc(reshape(norm_rf, [new_sizes(1), new_sizes(2)]))
                colormap(brewermap([],'RdBu'))
                caxis([0,1])
                axis square
                if nd == 1
                    title(num2str(datarun.cell_ids(cell_indices(rgc))))
                end
                               
                figure(12)
                subplot(3,1,nd)
                plot(V(:,nd))
                axis square
                axis tight
            end
        pause    
        end
    
        tmp_var_exp(rgc) = si(1);
    end
    SVDVarExpl{cts} = tmp_var_exp;
end


