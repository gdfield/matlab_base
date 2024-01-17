datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2013-10-28-0/data002/data002';
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-10-1-0.48-11111-60x60-60.35.xml';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

frame_num = 588;
[mov,height,width,duration,refresh] = get_movie(movie_path, datarun.triggers, frame_num);




