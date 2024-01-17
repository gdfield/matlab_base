%% set data paths

datapath_WT = get_cng_wn_datapaths('WT', 0);
datapath_1M = get_cng_wn_datapaths('1M', 0);
datapath_2M = get_cng_wn_datapaths('2M', 0);
datapath_3M = get_cng_wn_datapaths('3M', 0);
datapath_4M = get_cng_wn_datapaths('4M', 0);
datapath_5M = get_cng_wn_datapaths('5M', 0);
datapath_7M = get_cng_wn_datapaths('7M', 0);

%%
datapath = datapath_WT;
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

fraction_kept = zeros(length(datapath), 1);
svd_pass = zeros(length(datapath), 1);
var_pass = zeros(length(datapath), 1);
for dset = 1:length(datapath)
    [height, width, num_frames] = size(datarun{dset}.stas.stas{1});

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
    
    svd_pass(dset) = length(svd_indices) / length(datarun{dset}.cell_ids);
    var_pass(dset) = length(var_indices) / length(datarun{dset}.cell_ids);
    fraction_kept(dset) = length(keeper_indices) / length(datarun{dset}.cell_ids);
    
end

%% 
%datarun{dset} = load_sta(datarun{dset}, 'load_sta', datarun{1}.cell_ids(1:4));
%datarun{dset} = get_sta_summaries(datarun{dset}, datarun{1}.cell_ids(1:4));

% time course parameters to vary


clear init_params
init_params.fit_scale_one = true;
init_params.fit_scale_two = true;
init_params.fit_tau_one = true;
init_params.fit_tau_two = true;
init_params.initial_n_one_filters = 9;
init_params.initial_n_two_filters = 9;
init_params.fit_n_one_filters = false;
init_params.fit_n_two_filters = false;
init_params.verbose = false;
verbose = false;


for dset = 1:length(datapath)

    keeper_tcs = datarun{dset}.stas.keeper_tcs;
    keeper_indices = datarun{dset}.stas.keeper_indices;    
    sta_fits{length(datarun{dset}.cell_ids)} = [];
    datarun{dset}.matlab.sta_fits = sta_fits;

    clear tc_fit
    %fit_params = zeros(size(keeper_tcs,1), 7);
    for rgc = 1:size(keeper_tcs,1)
        [temp_fit, temp_params] = fit_time_course(keeper_tcs(rgc,:)', init_params);

        if verbose
            figure(1); clf;
            plot(keeper_tcs(rgc,:) ./ norm(keeper_tcs(rgc,:)), 'k')
            hold on
            plot(temp_fit ./ norm(temp_fit), 'r')
            hold off
            %temp_params
            %pause
        end
                
        tc_fit.scale_one = temp_params(1);
        tc_fit.scale_two = temp_params(2);
        tc_fit.tau_one = temp_params(3);
    	tc_fit.tau_two = temp_params(4);
        tc_fit.n_one_filters = temp_params(5);
        tc_fit.n_two_filters = temp_params(6);
        
        datarun{dset}.matlab.sta_fits{keeper_indices(rgc)} = tc_fit;
        

    end
end


time_to_zero = [];
time_to_peak = [];
cells_per_rec = zeros(length(datapath), 1);
for dset = 1:length(datapath)
    
    keeper_indices = datarun{dset}.stas.keeper_indices;
    [ttz, ttp] = get_time_to_zero(datarun{dset},...
                    datarun{dset}.cell_ids(keeper_indices),...
                    'step_size', 0.01, 'start_time', 2.5, 'end_time', 9.5,...
                    'ms_per_frame', 2/60.35);   
                
    time_to_zero = [time_to_zero, ttz];
    time_to_peak = [time_to_peak, ttp];
    cells_per_rec(dset) = length(keeper_indices);
end
frac_monophasic_tcs = length(find(isnan(time_to_zero) == 1)) ./ length(time_to_zero);


histogram(time_to_zero, 30);
histogram(time_to_peak, 30);

figure(1); clf;
hold on
clear inds_end inds_begin
inds_end = 0;
for dset = 1:length(cells_per_rec)
    inds_begin = 1 + inds_end;
    inds_end = sum(cells_per_rec(1:dset));
    [f, x] = ecdf(time_to_zero(inds_begin:inds_end));
    plot(x, f)
    pause
end
axis square
title('4M cdfs')
print(1, '~/Desktop/4M-cdfs.pdf', '-dpdf')


%% save WT
cd ~/Desktop/cng-tc-results
cng_tc_WT.time_to_zero = time_to_zero;
cng_tc_WT.time_to_peak = time_to_peak;
cng_tc_WT.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_WT.svd_pass = svd_pass;
cng_tc_WT.var_pass = var_pass;
cng_tc_WT.fraction_kept = fraction_kept;
cng_tc_WT.cells_per_rec = cells_per_rec;

save cng_tc_WT cng_tc_WT

%% save 1M
cd ~/Desktop/cng-tc-results
cng_tc_1M.time_to_zero = time_to_zero;
cng_tc_1M.time_to_peak = time_to_peak;
cng_tc_1M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_1M.svd_pass = svd_pass;
cng_tc_1M.var_pass = var_pass;
cng_tc_1M.fraction_kept = fraction_kept;
cng_tc_1M.cells_per_rec = cells_per_rec;

save cng_tc_1M cng_tc_1M

%% save 2M
cd ~/Desktop/cng-tc-results
cng_tc_2M.time_to_zero = time_to_zero;
cng_tc_2M.time_to_peak = time_to_peak;
cng_tc_2M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_2M.svd_pass = svd_pass;
cng_tc_2M.var_pass = var_pass;
cng_tc_2M.fraction_kept = fraction_kept;
cng_tc_2M.cells_per_rec = cells_per_rec;

save cng_tc_2M cng_tc_2M

%% save 3M
cd ~/Desktop/cng-tc-results
cng_tc_3M.time_to_zero = time_to_zero;
cng_tc_3M.time_to_peak = time_to_peak;
cng_tc_3M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_3M.svd_pass = svd_pass;
cng_tc_3M.var_pass = var_pass;
cng_tc_3M.fraction_kept = fraction_kept;
cng_tc_3M.cells_per_rec = cells_per_rec;

save cng_tc_3M cng_tc_3M
%% save 4M
cd ~/Desktop/cng-tc-results
cng_tc_4M.time_to_zero = time_to_zero;
cng_tc_4M.time_to_peak = time_to_peak;
cng_tc_4M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_4M.svd_pass = svd_pass;
cng_tc_4M.var_pass = var_pass;
cng_tc_4M.fraction_kept = fraction_kept;
cng_tc_4M.cells_per_rec = cells_per_rec;

save cng_tc_4M cng_tc_4M

%% save 5M
cd ~/Desktop/cng-tc-results
cng_tc_5M.time_to_zero = time_to_zero;
cng_tc_5M.time_to_peak = time_to_peak;
cng_tc_5M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_5M.svd_pass = svd_pass;
cng_tc_5M.var_pass = var_pass;
cng_tc_5M.fraction_kept = fraction_kept;
cng_tc_5M.cells_per_rec = cells_per_rec;

save cng_tc_5M cng_tc_5M


%% save 7M
cd ~/Desktop/cng-tc-results
cng_tc_7M.time_to_zero = time_to_zero;
cng_tc_7M.time_to_peak = time_to_peak;
cng_tc_7M.frac_monophasic_tcs = frac_monophasic_tcs;

cng_tc_7M.svd_pass = svd_pass;
cng_tc_7M.var_pass = var_pass;
cng_tc_7M.fraction_kept = fraction_kept;
cng_tc_7M.cells_per_rec = cells_per_rec;

save cng_tc_7M cng_tc_7M

%%
cd ~/Desktop/cng-tc-results
load cng_tc_WT
load cng_tc_1M
load cng_tc_2M
load cng_tc_3M
load cng_tc_4M
load cng_tc_5M
load cng_tc_7M

% plot fraction of cells passing the svd threshold
figure(1); clf;
plot(1, cng_tc_WT.svd_pass, 'kx')
hold on
errorbar(1, mean(cng_tc_WT.svd_pass), std(cng_tc_WT.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(2, cng_tc_1M.svd_pass, 'kx')
errorbar(2, mean(cng_tc_1M.svd_pass), std(cng_tc_1M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(3, cng_tc_2M.svd_pass, 'kx')
errorbar(3, mean(cng_tc_2M.svd_pass), std(cng_tc_2M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(4, cng_tc_3M.svd_pass, 'kx')
errorbar(4, mean(cng_tc_3M.svd_pass), std(cng_tc_3M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(5, cng_tc_4M.svd_pass, 'kx')
errorbar(5, mean(cng_tc_4M.svd_pass), std(cng_tc_4M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(6, cng_tc_5M.svd_pass, 'kx')
errorbar(6, mean(cng_tc_5M.svd_pass), std(cng_tc_5M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(7, cng_tc_7M.svd_pass, 'kx')
errorbar(7, mean(cng_tc_7M.svd_pass), std(cng_tc_7M.svd_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold off
axis([1 7 0 1])
axis square
title('SVD threshold')

% plot fraction of cells passing the variance threshold
figure(2); clf;
plot(1, cng_tc_WT.var_pass, 'kx')
hold on
errorbar(1, mean(cng_tc_WT.var_pass), std(cng_tc_WT.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(2, cng_tc_1M.svd_pass, 'kx')
errorbar(2, mean(cng_tc_1M.var_pass), std(cng_tc_1M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(3, cng_tc_2M.svd_pass, 'kx')
errorbar(3, mean(cng_tc_2M.var_pass), std(cng_tc_2M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(4, cng_tc_3M.svd_pass, 'kx')
errorbar(4, mean(cng_tc_3M.var_pass), std(cng_tc_3M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(5, cng_tc_4M.svd_pass, 'kx')
errorbar(5, mean(cng_tc_4M.var_pass), std(cng_tc_4M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(6, cng_tc_5M.svd_pass, 'kx')
errorbar(6, mean(cng_tc_5M.var_pass), std(cng_tc_5M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(7, cng_tc_7M.svd_pass, 'kx')
errorbar(7, mean(cng_tc_7M.var_pass), std(cng_tc_7M.var_pass), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold off
axis([1 7 0 1])
axis square
title('variance threshold')


% plot fraction of cells passing the both thresholds
figure(3); clf;
plot(1, cng_tc_WT.fraction_kept, 'kx')
hold on
errorbar(1, mean(cng_tc_WT.fraction_kept), std(cng_tc_WT.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(2, cng_tc_1M.svd_pass, 'kx')
errorbar(2, mean(cng_tc_1M.fraction_kept), std(cng_tc_1M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(3, cng_tc_2M.svd_pass, 'kx')
errorbar(3, mean(cng_tc_2M.fraction_kept), std(cng_tc_2M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(4, cng_tc_3M.svd_pass, 'kx')
errorbar(4, mean(cng_tc_3M.fraction_kept), std(cng_tc_3M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(5, cng_tc_4M.svd_pass, 'kx')
errorbar(5, mean(cng_tc_4M.fraction_kept), std(cng_tc_4M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(6, cng_tc_5M.svd_pass, 'kx')
errorbar(6, mean(cng_tc_5M.fraction_kept), std(cng_tc_5M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
plot(7, cng_tc_7M.svd_pass, 'kx')
errorbar(7, mean(cng_tc_7M.fraction_kept), std(cng_tc_7M.fraction_kept), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold off
axis([1 7 0 1])
axis square
title('both threshold')


%%

figure(4); clf;
[f_1M, x_1M] = ecdf(cng_tc_1M.time_to_zero);
plot(x_1M, f_1M, 'Color', '[0 0 1]')
hold on
pause
[f_2M, x_2M] = ecdf(cng_tc_2M.time_to_zero);
plot(x_2M, f_2M, 'Color', '[0.1 0.1 1]')
pause
[f_3M, x_3M] = ecdf(cng_tc_3M.time_to_zero);
plot(x_3M, f_3M, 'Color', '[0.2 0.2 1]')
pause
[f_4M, x_4M] = ecdf(cng_tc_4M.time_to_zero);
plot(x_4M, f_4M, 'Color', '[0.5 0.5 1]')
pause
[f_5M, x_5M] = ecdf(cng_tc_5M.time_to_zero);
plot(x_5M, f_5M, 'Color', '[0.6 0.6 1]')
pause
[f_7M, x_7M] = ecdf(cng_tc_7M.time_to_zero);
plot(x_7M, f_7M, 'Color', '[0.7 0.7 1]')
pause
[f_WT, x_WT] = ecdf(cng_tc_WT.time_to_zero);
plot(x_WT, f_WT, 'r')
pause
hold off












