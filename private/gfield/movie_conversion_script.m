cd ~/Development/natural-movies
movie_path = '~/Development/natural-movies/original/';
save_path = '~/Development/natural-movies/conversion/';
movie_name = 'trellis_mr';
movie_ext = '.mp4';
path_name = [movie_path movie_name movie_ext];
reduction_factor = 1; % fraction of movie to use


vidObj = VideoReader(path_name);

desired_width = 600;
desired_height = 600;

num_frames = floor(vidObj.Duration * vidObj.FrameRate);
temp_height = vidObj.Height;
temp_width = vidObj.Width;

num_frames_to_write = floor(num_frames * reduction_factor) -1;


height_begin = floor((temp_height - desired_height) / 2);
height_end = temp_height - height_begin - 1;
width_begin = floor((temp_width - desired_height) / 2);
width_end = temp_width - width_begin -1;

mov = zeros(desired_height, desired_width, num_frames_to_write);


% height_begin = temp_height - desired_height+1;
% height_end = temp_height;
% width_begin = 1;
% width_end = desired_width;


for fm = 1:num_frames_to_write
    tmp_frame = read(vidObj,fm);
    tmp_frame = double(tmp_frame);
    tmp_fm = sum(tmp_frame, 3) ./ 3;
    mov(:,:,fm) = tmp_fm(height_begin:height_end,width_begin:width_end);
    
    show_image = repmat(tmp_fm(height_begin:height_end,width_begin:width_end), [1 1 3]);
    image(uint8(show_image))
    drawnow
end
mov = uint8(mov);
size(mov)

save([save_path,movie_name], 'mov', '-v7.3')

write_movie([save_path, movie_name, '.mat'], [save_path, movie_name, '.rawMovie'], 1);


%% make dark version

for fm = 1:num_frames_to_write
    tmp_frame = read(vidObj,fm);
    tmp_frame = double(tmp_frame);
    tmp_fm = sum(tmp_frame, 3) ./ 3;
    
    tmp_fm = tmp_fm ./ 3;
    
    mov(:,:,fm) = tmp_fm(height_begin:height_end,width_begin:width_end);
    
    show_image = repmat(tmp_fm(height_begin:height_end,width_begin:width_end), [1 1 3]);
    image(uint8(show_image))
    drawnow
end

dark_name = [movie_name,'_dark'];
save([save_path,dark_name], 'mov', '-v7.3')
write_movie([save_path, dark_name, '.mat'], [save_path, dark_name, '.rawMovie'], 1);

%% make bright version

for fm = 1:num_frames_to_write
    tmp_frame = read(vidObj,fm);
    tmp_frame = double(tmp_frame);
    tmp_fm = sum(tmp_frame, 3) ./ 3;
    
    tmp_fm = tmp_fm ./ 3 + 169;
    
    mov(:,:,fm) = tmp_fm(height_begin:height_end,width_begin:width_end);
    
    show_image = repmat(tmp_fm(height_begin:height_end,width_begin:width_end), [1 1 3]);
    image(uint8(show_image))
    drawnow
end

bright_name = [movie_name,'_bright'];
save([save_path,bright_name], 'mov', '-v7.3')
write_movie([save_path, bright_name, '.mat'], [save_path, bright_name, '.rawMovie'], 1);


%% make alternating version

lum_offset = zeros(desired_height, desired_width, num_frames_to_write);
lum_offset(:,:,121:240) = 169;
lum_offset(:,:,361:406) = 169;

for fm = 1:num_frames_to_write
    tmp_frame = read(vidObj,fm);
    tmp_frame = double(tmp_frame);
    tmp_fm = sum(tmp_frame, 3) ./ 3;
    
    tmp_fm = tmp_fm ./ 3;
    
    mov(:,:,fm) = tmp_fm(height_begin:height_end,width_begin:width_end) + lum_offset(:,:,fm);
    
    show_image = repmat(tmp_fm(height_begin:height_end,width_begin:width_end) + lum_offset(:,:,fm), [1 1 3]);
    image(uint8(show_image))
    drawnow
end

alt_name = [movie_name,'_alt'];
save([save_path,alt_name], 'mov', '-v7.3')
write_movie([save_path, alt_name, '.mat'], [save_path, alt_name, '.rawMovie'], 1);

