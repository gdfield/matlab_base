%%%%%%%% PHOTOGRAPHIC MAPPING %%%%%%%%%


% choose dataset
datarun = load_data('2009-04-13-5/data005/data005');


% load datarun
datarun = load_index(datarun, 'basepath', '/Volumes/Archive/');
datarun = load_sta(datarun,'load_sta',[]);
datarun = load_params(datarun);
datarun = load_ei(datarun,[]);


% load alignment information from disk
datarun = load_monitor_alignment(datarun);


% click alignment points
datarun = compute_monitor_to_array_transformation(datarun);


% verify the alignment points look good
verify_array_monitor_alignment(datarun)
figure;plot_rf_summaries(datarun,'all','label',true,'array',true,'array_label',true)


% save the result
save_monitor_alignment(datarun);
