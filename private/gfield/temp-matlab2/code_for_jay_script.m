%% code to load in movie

datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-06-24-0/data005_MR/data005_MR';
moviepath = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';

datarun = load_data(datapath);
datarun = load_neurons(datarun);

num_frames = 50000;
[mov,height,width,duration,refresh] = get_movie(moviepath, datarun.triggers, num_frames);
mov = squeeze(mov(:,:,1,:)) - 0.5;

% mov is now X x Y x frames and varies from -0.48 to + 0.48

%% Get SNLs, fits and plot

datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/...';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');

marks_params.thresh = 4.5;
cell_types = {'ON type1', 'OFF type1','ON type2', 'OFF type2', 'ON type3', 'OFF type3'};
num_types = length(cell_types);
datarun = get_sta_summaries(datarun, cell_types, 'marks_params', marks_params);

datarun = load_java_movie(datarun, 'path_to_movie');

datarun = get_snls(datarun, {'cell_type1', 'cell_type2'}, 'stimuli', 400000,'verbose', true, 'new', true, 'fit', 'cum_norm');

type1_indices = get_cell_indices(datarun, {'cell_type1'});
figure(3); clf;
hold on

% plot the snls fits
for rgc = 1:length(type1_indices)

    fit_params = datarun.stas.snls{type1_indices(rgc)}.fit_params;
    fitx = linspace(start_plot, end_plot, num_points);
    fit_fcn = @(x)(fit_params.a * normcdf(fit_params.b * x - fit_params.c, 0, 1));
    plot(fitx,fit_fcn(fitx),'b-')
 
end   
