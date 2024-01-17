%% Local vs. distant return
visionDataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 
datarun = load_data(visionDataPath);  
datarun = load_params(datarun); 
datarun = load_sta(datarun); 
datarun = load_ei(datarun,'all'); 
cell_ids = [125 137 156 198 302 542 544 557 736 901 1085 1159 1323 1516 ...
    1657 1862 1881 2869 4097 4188 4217 4247 4399 4414 4757 4774 4816 ...
    4894 4968 5152 5238 5326 5371 5388 5672 5748 6001 6139 6256 6408 ...
    6678 6739 6978 6995 7144 7202];

plot_rf_summaries(datarun, cell_ids, 'clear', false, 'label', false,...
        'plot_fits', true, 'fit_color', [0.8 0.8 0.8],...
        'fill_color',[0.8 0.8 0.8]);
axis off;  title('ON Ps analyzed for electrical stimulation response'); 

% For each neuron, flip through and determine which stimulation electrodes were analyzed from elecResp files
