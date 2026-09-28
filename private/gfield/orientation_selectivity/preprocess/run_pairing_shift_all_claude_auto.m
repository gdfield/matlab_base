%% run_pairing_shift_all_claude_auto.m
% Runs os_pairing_shift_test_claude_auto on every dataset in new_os_datasets()
% and writes a results table. Set ps_subset = {'<dataset>', ...} in the
% workspace first to run only some datasets; results are merged by dataset
% name into the existing results file.
% 2026-09-27 GDF + Claude

out_root = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
res_mat  = [out_root, 'pairing_shift_test_claude_auto.mat'];
res_csv  = [out_root, 'pairing_shift_test_claude_auto.csv'];
dlist = new_os_datasets();
if exist(res_mat, 'file'), load(res_mat, 'PS'); else, PS = struct([]); end
for i = 1:numel(dlist)
    dn = regexp(dlist(i).grating_datapath, '\d{4}-\d{2}-\d{2}-\d', 'match', 'once');
    if exist('ps_subset', 'var') && ~isempty(ps_subset) && ~ismember(dn, ps_subset), continue; end
    fprintf('pairing-shift test: %s\n', dn);
    R = os_pairing_shift_test_claude_auto(dlist(i), 'n_perm', 1000);
    R.run_date = datestr(now, 'yyyy-mm-dd HH:MM');
    if isempty(PS), PS = R; else
        j = find(strcmp({PS.dataset}, dn));
        if isempty(j), PS(end+1) = R; else, PS(j) = R; end
    end
    save(res_mat, 'PS');
end
T = table({PS.dataset}', [PS.n_trials]', [PS.n_reps]', [PS.n_cells]', [PS.observed]', ...
    [PS.max_other_shift]', [PS.best_shift]', [PS.null_mean]', [PS.null_p95]', [PS.null_max]', [PS.p_perm]', ...
    'VariableNames', {'dataset','n_trials','n_reps','n_cells','observed','max_other_shift', ...
    'best_shift','null_mean','null_p95','null_max','p_perm'});
writetable(T, res_csv);
