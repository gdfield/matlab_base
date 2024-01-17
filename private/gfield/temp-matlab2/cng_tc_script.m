%% set data paths

datapath_WT = get_cng_wn_datapaths('WT', 0);
datapath_1M = get_cng_wn_datapaths('1M', 0);
datapath_2M = get_cng_wn_datapaths('2M', 0);
datapath_3M = get_cng_wn_datapaths('3M', 0);
datapath_4M = get_cng_wn_datapaths('4M', 0);
datapath_5M = get_cng_wn_datapaths('5M', 0);
datapath_7M = get_cng_wn_datapaths('7M', 0);

%%
datapath = datapath_4M;
% load data and get STA summaries
for dset = 1:length(datapath)
    datarun{dset} = load_data(datapath{dset});
    datarun{dset} = load_neurons(datarun{dset});
    datarun{dset} = load_params(datarun{dset});
    datarun{dset} = load_sta(datarun{dset}, 'load_sta', 'all');
    marks_params.thresh = 4.0;
    datarun{dset} = get_sta_summaries(datarun{dset}, 'all', 'marks_params', marks_params);
    num_rgcs(dset) = length(datarun{dset}.cell_ids);


    % set filter radius for reducing STA noise for SVD
    filt_params.radius = 0.75;

    % decide if you want to plot shit to check out SVD performance
    verbose = 0; % if 1, then plot

    % get STA size and initialize several variable to speed code.
    [height, width, num_frames] = size(datarun{dset}.stas.stas{1});
    svd_tc = zeros(num_rgcs(dset), num_frames);
    rank_one_var = zeros(num_rgcs(dset),1);
    svd_rf = zeros(num_rgcs(dset), height*width);
    
    % loop over RGCs
    for rgc = 1:num_rgcs(dset)

        % filter the STA
        temp_sta = squeeze(datarun{dset}.stas.stas{rgc});

        
        % get center of mass and of STA and cut out an ROI to reduct dim of stim
        if datarun{dset}.vision.sta_fits{rgc}.mean == [-1 -1] % handle case where there isn't a good fit in Vision
            temp_com = [round(height/2), round(width/2)];
        else
            temp_com = abs([0 40] - round(datarun{dset}.vision.sta_fits{rgc}.mean)) +[1 1];
        end

        if temp_com(2) > height
            temp_com(2) = round(height./2);
        end
        if temp_com(1) > width
            temp_com(1) = round(width./2);
        end

        if isnan(datarun{dset}.vision.sta_fits{rgc}.sd) % handle case where there isn't a good fit in Vision
            temp_size = 4;
        else
            temp_size = 4*sqrt(prod(datarun{dset}.vision.sta_fits{rgc}.sd));
        end

        window_size = ceil(temp_size);
        bg_height = floor(temp_com(2) - window_size);
        end_height = floor(temp_com(2) + window_size);
        bg_width = floor(temp_com(1) - window_size);
        end_width = floor(temp_com(1) + window_size);
        if bg_height < 1; bg_height = 1; end
        if end_height > height; end_height = height; end
        if bg_width < 1; bg_width = 1; end
        if end_width > width; end_width = width; end
        windowed_height = end_height - bg_height +1;
        windowed_width = end_width - bg_width +1 ;
        temp_sta = temp_sta(bg_height:end_height, bg_width:end_width,:);
        
        filt_sta = zeros(size(temp_sta));
        for fm = 1:num_frames
            filt_fm = rf_filtered(squeeze(temp_sta(:,:,fm)), filt_params);
            filt_sta(:,:,fm) = filt_fm;
        end

        % do SVD
        reshape_sta = reshape(filt_sta, [windowed_height * windowed_width, num_frames]);
        [U, S, V] = svd(reshape_sta);   

        % plot shit
        if verbose
            figure(1); clf;
            subplot(1,3,1)
            plot((diag(S).^2) ./sum(diag(S).^2), 'ok')
            axis square
            axis tight
            ylabel('var exp')
            xlabel('rank')
            subplot(1,3,2)
            image(norm_image(reshape(U(:,1), [windowed_height, windowed_width])))
            title(num2str(datarun{dset}.cell_ids(rgc)))
            axis square
            axis off
            subplot(1,3,3)
            plot(V(:,1), 'k')
            hold on
            plot([0,num_frames], [0 0], 'k')
            axis tight
            axis square
            xlabel('frames')
            drawnow
        end

        % store output of TC, RF and rank 1 variance explained.
        svd_tc(rgc,:) = V(:,1)'; 
        %svd_rf(rgc,:) = U(:,1)';
        rank_one_var(rgc) = (diag(S(1)).^2) ./sum(diag(S).^2);

    end

    datarun{dset}.stas.svd_tc = svd_tc;
    %datarun{dset}.stas.svd_rf = svd_rf;
    datarun{dset}.stas.rank_one_var = rank_one_var;
end

%%
% checkout a histogram of the rank-1 variance explained.

for dset = 1:length(datapath)
    figure(1); clf;
    hist(datarun{dset}.stas.rank_one_var, 50)
    pause
    
    figure(1); clf;
    hist(var(datarun{dset}.stas.svd_tc'))
    pause
    
    % set thresholds for variance explained and for variance in the TC
    % good for WT at NDF0
    svd_threshold = 0.6; 
    tc_var_threshold = 0.025; 
    
    % good for WT at NDF2
    %svd_threshold = 0.25;
    %tc_var_threshold = 0.05;

    % get indices that comply with both thresholds
    svd_indices = find(datarun{dset}.stas.rank_one_var > svd_threshold);
    var_indices = find(var(datarun{dset}.stas.svd_tc') > tc_var_threshold);
    keeper_indices = intersect(svd_indices, var_indices);

    % extract the TCs to keep
    keeper_tcs = datarun{dset}.stas.svd_tc(keeper_indices,:);
    offset_vals = mean([keeper_tcs(:,1), keeper_tcs(:,num_frames)]');
    keeper_tcs = keeper_tcs - repmat(offset_vals', [1,num_frames]);

    % sanity check: plot the TCs and make sure they look 'good'
    figure(1);clf;
    plot(keeper_tcs')
    title([num2str(size(keeper_tcs, 1)),' quality TCs'])
    pause
    
    datarun{dset}.stas.keeper_indices = keeper_indices;
    datarun{dset}.stas.keeper_tcs = keeper_tcs;
    
end



%% Get Power Spectra for all of the RGCs
samplerate = 60.35;
num_keeper_rgcs = length(keeper_indices);
rgc_ps = zeros(num_keeper_rgcs, num_frames/2 +1);
for rgc = 1:num_keeper_rgcs
    test_tc = keeper_tcs(rgc,:);
    test_tc = test_tc ./ norm(test_tc);
    [powerspec_xvalues, mean_powerspec, fft_signal] = PowerSpectrumFinder(test_tc,samplerate);
    freq_ranges = powerspec_xvalues;
    norm_ps = mean_powerspec/max(mean_powerspec); % normalize these values to the peak
    rgc_ps(rgc, :) = norm_ps;
end

% compute max TF tuning
[max_TF, max_inds] = max(rgc_ps, [], 2);
bp_attenuation = 1 - rgc_ps(:,1);
temp_hist = histogram2(powerspec_xvalues(max_inds), bp_attenuation', 10, 'DisplayStyle', 'bar3');
xlabel('frequencies')
ylabel('attenuation')


hist(max(rgc_ps, [], 2), 40)
semilogy(powerspec_xvalues, rgc_ps)

%% Inspect the PCs of the TCs

num_dims_to_plot = 12;
for dset = 1:length(datapath)
    keeper_tcs = datarun{dset}.stas.keeper_tcs;
    % normalize the tcs to have unit length
    for rgc = 1:size(keeper_tcs,1)
        modified_tcs(rgc,:) = keeper_tcs(rgc,:) ./ norm(keeper_tcs(rgc,:));
        %modified_tcs(rgc,:) = modified_tcs(rgc,:) ./ ext(modified_tcs(rgc,:));
    end

    [tc_pcs, pc_weights, tc_eigenvalues] = pca(modified_tcs);
    figure(dset);clf;
    subplot(1,2,1)
    semilogy(tc_eigenvalues(1:num_dims_to_plot).^2 ./ sum(tc_eigenvalues.^2), 'ko');
    xlabel('dimensions')
    ylabel('var explained')
    title('PCA on Time Courses')
    axis square
    subplot(1,2,2)
    plot(tc_pcs(:,1), 'k')
    hold on
    plot(tc_pcs(:,2), 'k--')
    plot(tc_pcs(:,3), 'k-.')
    plot(tc_pcs(:,4), 'k:')
    %plot(tc_pcs(:,5), 'r')
    plot([0, length(tc_pcs(:,1))], [0, 0], 'k')
    title('Time Course PCs')
    xlabel('frames')
    ylabel('amplitude')
    axis square
    axis tight
    hold off
    
    datarun{dset}.stas.tc_pcs = tc_pcs;
    datarun{dset}.stas.tc_eigenvalues = tc_eigenvalues;
    datarun{dset}.stas.pc_weights = pc_weights;
end



figure(5); clf;
for ndim = 1:4
    subplot(2,2,ndim)
    plot(datarun{1}.stas.tc_pcs(:,ndim), 'k');
    hold on
    plot(datarun{2}.stas.tc_pcs(:,ndim), 'g');
    plot(datarun{3}.stas.tc_pcs(:,ndim), 'c');
    plot(datarun{4}.stas.tc_pcs(:,ndim), 'b');
    hold off
    title(['PC ', num2str(ndim)])
end


cd ~/Desktop/
datarun_4M = datarun;
save datarun_4M datarun_4M



