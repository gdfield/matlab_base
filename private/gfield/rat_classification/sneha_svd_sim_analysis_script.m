%% SNR check on long STA
sta_long = sim_sta_offt4_cell3_306187frames;
% check SNR
noise_block = sta_long(31:40,66:75,:);
noise_est = std(noise_block(:));
sorted_sta = sort(abs(sta_long(:)),'descend');
signal_est = mean(sorted_sta(1:5));
long_snr = signal_est/noise_est

%% SNR check on short STA
sta_short = sim_sta_offt4_cell3_212238frames;
% check SNR
noise_block = sta_short(31:40,66:75,:);
noise_est = std(noise_block(:));
sorted_sta = sort(abs(sta_short(:)),'descend');
signal_est = mean(sorted_sta(1:5));
short_snr = signal_est/noise_est

%% Inspect the STAs

% long STA
temp_sta = sta_long;
[~, ext_ind] = max(abs(temp_sta(:)));
ext_val = temp_sta(ext_ind);
norm_sta = temp_sta ./ abs(ext_val);
norm_long_sta = squeeze((norm_sta + 1) ./ 2);

% short STA
temp_sta = sta_short;
[~, ext_ind] = max(abs(temp_sta(:)));
ext_val = temp_sta(ext_ind);
norm_sta = temp_sta ./ abs(ext_val);
norm_short_sta = squeeze((norm_sta + 1) ./ 2);


num_frames = size(norm_short_sta, 3);
for fm = 1:num_frames
    subplot(2,1,1)
    image(repmat(norm_long_sta(:,:,fm), [1,1,3]))
    subplot(2,1,2)
    image(repmat(norm_short_sta(:,:,fm), [1,1,3]))
    pause
end


%% perform SVD on simulated RF

frame_start = 14;
frame_end = 30;

model_sta = reshape(sta_long, [3200, 30]);
model_sta = model_sta(:,frame_start:frame_end);
[U, D, V] = svd(model_sta);
long_spectrum = diag(D);
long_si = long_spectrum.^2 ./ sum(long_spectrum.^2);

model_sta = reshape(sta_short, [3200, 30]);
model_sta = model_sta(:,frame_start:frame_end);
[U, D, V] = svd(model_sta);
short_spectrum = diag(D);
short_si = short_spectrum.^2 ./ sum(short_spectrum.^2);

figure(1); clf;
plot(long_si, 'ko')
hold on
plot(short_si, 'ro')
axis auto

figure(2); clf;
plot(cumsum(long_si), 'k')
hold on
plot(cumsum(short_si), 'r')



model_space_kernal = U(:,1);
model_space_kernal = reshape(model_space_kernal, [40, 80]);
% normalize STA range
[~, ext_ind] = max(abs(model_space_kernal(:)));
ext_val = model_space_kernal(ext_ind);
norm_smodel_space_kernal = model_space_kernal ./ abs(ext_val);

norm_smodel_space_kernal = (norm_smodel_space_kernal + 1) ./ 2;

image(repmat(norm_smodel_space_kernal, [1,1,3]))



%% inspect STA

cell_sta = STA_data_offt4_3;
[~, ext_ind] = max(abs(cell_sta(:)));
ext_val = cell_sta(ext_ind);
norm_cell_sta = cell_sta ./ abs(ext_val);
norm_cell_sta = squeeze((norm_cell_sta + 1) ./ 2);

num_frames = size(norm_cell_sta, 3);
for fm = 1:num_frames
    subplot(3,1,1)
    image(repmat(norm_long_sta(:,:,fm), [1,1,3]))
    axis off
    subplot(3,1,2)
    image(repmat(norm_short_sta(:,:,fm), [1,1,3]))
    axis off
    subplot(3,1,3)
    image(repmat(norm_cell_sta(:,:,fm), [1,1,3]))
    axis off
    pause
end

cell_sta = reshape(cell_sta, [3200, 30]);
cell_sta = cell_sta(:, frame_start:frame_end);
[U, D, V] = svd(cell_sta);
data_spectrum = diag(D);
data_si = data_spectrum.^2 ./ sum(data_spectrum.^2);

plot(data_si, 'bo')
plot(cumsum(data_si), 'b')


space_kernal = U(:,1);
space_kernal = reshape(space_kernal, [40, 80]);

% normalize STA range
[~, ext_ind] = max(abs(space_kernal(:)));
ext_val = space_kernal(ext_ind);
norm_space_kernal = space_kernal ./ abs(ext_val);

norm_space_kernal = (norm_space_kernal + 1) ./ 2;

image(repmat(norm_space_kernal, [1,1,3]))

% check snr on data
cell_sta = STA_data_offt4_3;
noise_block = cell_sta(31:40,66:75,1:10);
noise_est = robust_std(noise_block(:));
sorted_sta = sort(abs(cell_sta(:)),'descend');
signal_est = mean(sorted_sta(1:5));
snr = signal_est/noise_est




