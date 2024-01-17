% load data, compute RFs and TCs
%2012-10-10-1
datarun = load_data('/Users/gdfield/Analysis/2012-10-10-1/data001-3600-7200s/data001-map/data001-map');
offt3 = [512 1368 2750 4188 4487 5732 6170 6511 6588];
offt5 = [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353,3874];

%2012-10-15-0
datarun = load_data('/Users/gdfield/Analysis/2012-10-15-0/data000-map/data000-map');
offt3 = [62,991,1156,4234,4278,4487,5733,6286,6931];
offt5 = [1246,2253,3695,5116,6260];

%2012-10-31
datarun = load_data('/Users/gdfield/Analysis/2012-10-31-0/data000-map/data000-map');
offt3 = [1577,1954,2659,7324,6033,4383,2884];
offt5 = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376];


datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
marks_params.thresh = 4.5;
datarun = get_sta_summaries(datarun, 'all','marks_params', marks_params);

% book keeping
%cell_type_num = get_cell_type_nums(datarun, 'on type1');
cell_indices = get_cell_indices(datarun, offt5);

num_rgcs = length(cell_indices);
for rgc = 1:num_rgcs

    % get the space-time STA
    temp_sta = datarun.stas.stas{cell_indices(rgc)};

    % squeeze and reshape the sta for SVD
    temp_sta = squeeze(temp_sta);
    temp_sta = reshape(temp_sta, [3200,30]);

    % perform SVD
    [U, S, V] = svd(temp_sta);
    
    [W, H] = nnmf(temp_sta,4);

    % plot stuff
    num_dims = 4;
    spaces = zeros(40, 80, num_dims);
    temporals = zeros(30, num_dims);
    score_thresh = 3.5;
    for dims = 1:num_dims
        spaces(:, :, dims) = reshape(U(:,dims), [40, 80]);
        figure(1)
        subplot(2, 2, dims)
        colormap jet
        imagesc(spaces(:,:, dims))

        temporals(:,dims) = V(:,dims);
        figure(2)
        subplot(2,2,dims)
        plot(temporals(:,dims))
        
        figure(3)
        plot(diag(S), 'o')
        
        
        figure(4)       
        S = diag(S);
        hold on
        plot(S(1)./S(15), S(2)./S(15), 'bo')
        hold off
        
        %compute robust std
        temp_r_sd = robust_std(U(:,dims));
        % put spatial component in z score
        zU = U(:,dims) ./ temp_r_sd;
        % threshold by z_scores
        zU(abs(zU)<score_thresh) = 0;
        
%         figure(4)
%         subplot(3,3,dims)
%         thresh_rf(:,:,dims) = reshape(abs(zU), [40, 80]);
%         imagesc(thresh_rf(:,:,dims)); colormap gray;
        
        
    end
    
%     figure(6)
%     size(thresh_rf)
%     cumulative_rf(:,:,dims) = sum(thresh_rf(:,:,1:3), 3);
%     imagesc(squeeze(cumulative_rf(:,:,dims))); colormap gray

    drawnow
    pause
end

