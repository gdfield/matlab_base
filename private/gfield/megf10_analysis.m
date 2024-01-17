I_cells = MIJ.getCurrentImage;

% switch image
I_processes = MIJ.getCurrentImage;

% switch image
I_control = MIJ.getCurrentImage;


I_cells = double(I_cells);
I_processes = double(I_processes);
I_control = double(I_control);


I_cells = I_cells - mean(I_cells(:));
I_processes = I_processes - mean(I_processes(:));


% take cross-cor
within_corr = xcorr2(I_cells, I_processes);
figure(1);
imagesc(within_corr)

%% recal with flipped image

rev_indices = [1024:-1:1];
auto_rotated = xcorr2(I_cells, I_processes(rev_indices,rev_indices));
figure(2)
imagesc(auto_rotated)

%% subtract the two

figure(3)
diff_image = within_corr - auto_rotated;
imagesc(diff_image)

figure(4)
plot(diff_image(1024,:)./max(auto_rotated(:,1024)) )

%%
% take cross-cor
control_corr = xcorr2(I_cells, I_control);
figure(11);
imagesc(control_corr)

figure(12)
control_diff = within_corr - control_corr;
imagesc(control_diff)

figure(13)
plot(control_diff(1024,:)./max(within_corr(:,1024)) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%filter parameters
Gauss_params.x_dim = 100;
Gauss_params.y_dim = 100;
Gauss_params.sd_x = 15;
Gauss_params.sd_y = 15;
Gauss_params.center_point_x = 50;
Gauss_params.center_point_y = 50;
temp_Gauss = make_Gaussian_two_d(Gauss_params);
temp_Gauss = temp_Gauss ./ max(temp_Gauss(:));

Gauss_area = Gauss_params.sd_x^2 *2*pi ./ 1e6;
I_cells_filt = conv2(I_cells, temp_Gauss, 'same') ./ Gauss_area;
I_processes_filt = conv2(I_processes, temp_Gauss, 'same') ./ Gauss_area;

figure(11)
imagesc(I_cells_filt);
figure(12)
imagesc(I_processes_filt);




% take cross-cor
within_corr = xcorr2(I_cells_filt*-1, I_processes);
figure(11);
imagesc(within_corr)

%% recal with flipped image

rev_indices = [1024:-1:1];
auto_rotated = xcorr2(I_cells_filt *-1, I_processes(rev_indices,rev_indices));
figure(12)
imagesc(auto_rotated)

%% subtract the two

figure(13)
diff_image = within_corr - auto_rotated;
imagesc(diff_image)

figure(14)
plot(diff_image(1024,:)./max(auto_rotated(:,1024)) )






