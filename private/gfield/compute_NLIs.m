function [NLIs, NLIs_err, NLIs_alt, NLIs_alt_err] = compute_NLIs(datarun, cell_type)

snl_range = 1.5;
num_points = 100;
num_datasets = length(datarun);


for dset = 1:num_datasets
    tmp_indices = get_cell_indices(datarun{dset}, cell_type);
    tmp_NLIs = zeros(length(tmp_indices), 1);
    tmp_NLIs_alt = zeros(length(tmp_indices), 1);  %
    
    for rgc = 1:length(tmp_indices)
        fit_params = datarun{dset}.stas.snls{tmp_indices(rgc)}.fit_params;
        fitx = linspace(-1*snl_range, snl_range, num_points);
        fitfcn = fit_params.a * normcdf(fit_params.b * fitx - fit_params.c, 0, 1);
        
        slope_at_zero = abs(diff(fitfcn(50:51)));
        slope_at_max = max(diff(fitfcn));

        tmp_NLIs(rgc) = log(slope_at_max ./ slope_at_zero);      
        tmp_NLIs_alt(rgc) = log(sum(fitfcn(51:100))/sum(fitfcn(1:50))); %
        
    end   
    tmp_NLIs = tmp_NLIs(~isnan(tmp_NLIs));
    tmp_NLIs = tmp_NLIs(isfinite(tmp_NLIs));
    NLIs(dset) = mean(tmp_NLIs);
    NLIs_err(dset) = std(real(tmp_NLIs)) ./ sqrt(length(tmp_NLIs));
    
    tmp_NLIs_alt = tmp_NLIs_alt(~isnan(tmp_NLIs_alt));  %
    tmp_NLIs_alt = tmp_NLIs_alt(isfinite(tmp_NLIs_alt));  %
    NLIs_alt(dset) = mean(tmp_NLIs_alt);  %
    NLIs_alt_err(dset) = std(real(tmp_NLIs_alt)) ./ sqrt(length(tmp_NLIs_alt)); %    
end