temp_spec.x_start = 0;
temp_spec.y_start = 0;
temp_spec.x_end = 800;
temp_spec.y_end = 600;
temp_spec.orientation = 0;
grating_duration = 8;
frame_rate = 60;

temp_spec.spatial_period = 60;
temp_spec.temporal_period = 60;
temp_spec.spatial_phase = 0;

[temp_frame, temp_tscale] = calc_reversing_grating_frame_intensities(temp_spec);
temp_frame = repmat(temp_frame, [1 1 3]);
temp_frame = temp_frame./2 + 0.5;

stimulus_set.max_frame = temp_frame;
stimulus_set.tscale = repmat(temp_tscale, 1, grating_duration * frame_rate / temp_spec.temporal_period);

h = figure('MenuBar','none', 'Color', [1 1 1], 'Position',[100 100 500 500]);

clear mov
num_frames = 240;
for fm = 1:num_frames
    
    clf

    temp_spec.spatial_period = 120;
    temp_spec.temporal_period = 120;
    temp_spec.spatial_phase = fm*2;

    [temp_frame, temp_tscale] = calc_reversing_grating_frame_intensities(temp_spec);
    temp_frame = repmat(temp_frame, [1 1 3]);
    temp_frame = temp_frame./2 + 0.5;

    
    image(temp_frame);
    axis equal
    axis off

    mov(fm) = getframe(h);
    
end

save_path = '~/Desktop/sta-movie';
movie2avi(mov,save_path,'FPS',60);