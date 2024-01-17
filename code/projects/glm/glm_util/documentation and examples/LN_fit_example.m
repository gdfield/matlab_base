%%{
clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
dsave_WN = '/Volumes/Lab/Users/Nora/LNFits/RPE/201602171/WN/OffPar/'; mkdir(dsave_WN);
dsave_NSEM = '/Volumes/Lab/Users/Nora/LNFits/RPE/201602171/NSEM/OffPar/';mkdir(dsave_NSEM);
cells = get_cell_ids(datarun_class, 'Off Parasol'); % cell ids to fit
datarun_class = load_sta(datarun_class,'load_sta', cells);
%%
%{
test_data = 'data027';
fit_data = 'data026';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1200, 30, 'cell_spec', cells,'visual_check', 1, 'stimulus_name', 'BW-8-1', 'seed', 22222);
FG = ln_fit_from_WN(cells, [Analysis_Path '/' fit_data '/' fit_data], 'BW-8-1-0.48-11111', 'testmovie', repeats.testmovie, 'testspikes', repeats.testspikes, 'stim_length', 1800, 'd_save', dsave_WN);
disp('done with WN glm fit')
%}
%% NSEM
%%{
test_data = 'data022';
fit_data = 'data025';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1200, 30,'cell_spec', cells,'visual_check', 0);
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
monitor_refresh = 100/median(diff(fit_datarun.triggers));
load('/Volumes/Lab/Users/Nora/Stimulus/NSEM_Movies/downsampledNSinterval.mat')
testmovie = fitmovie(:,:,1:1200);

%%
for i = 1
    disp(i)
    glm_cellinfo.cid           = cells(i);
    glm_cellinfo.cell_savename = num2str(cells(i));
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100);
    center(1) = 40 -datarun_class.vision.sta_fits{master_idx}.mean(2);
    center(2) = datarun_class.vision.sta_fits{master_idx}.mean(1);
    center = ceil(center);
    STA = squeeze(sum(datarun_class.stas.stas{master_idx},3));
    fittedLN     = ln_fit(fitspikes, fitmovie, center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA, 'iterations', 1);
    fittedLN = ln_predict(fittedLN, testmovie,'testspikes', repeats.testspikes(:,i));
    close all
    plotfilters_LN(fittedLN)
    plotrasters(fittedLN.xval, fittedLN.orig_glm);
end

%}
