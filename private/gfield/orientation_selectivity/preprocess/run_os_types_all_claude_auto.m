function done = run_os_types_all_claude_auto(budget_s)
% run_os_types_all_claude_auto  Resumable batch of os_types_hsimple_claude_auto over
% every dataset in os_cells_master_table.mat; skips datasets whose types.mat is
% version >= 2 (or that have v2_ERROR.txt); stops starting new datasets after budget_s seconds.
% Appends one line per dataset to os_analysis/claude_types/run_log.txt.
% 2026-09-27 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); dsl = unique(S.OS.dataset, 'stable');
if ~exist([oa 'claude_types'], 'dir'), mkdir([oa 'claude_types']); end
lg = fopen([oa 'claude_types/run_log.txt'], 'a'); t0 = tic; done = 0;
for k = 1:numel(dsl)
    tf = [oa 'claude_types/' dsl{k} '/types.mat'];
    if exist(tf, 'file'), q = load(tf, 'info'); if isfield(q.info, 'version') && q.info.version >= 2, done = done + 1; continue; end, end
    if exist([oa 'claude_types/' dsl{k} '/v2_ERROR.txt'], 'file'), done = done + 1; continue; end
    if toc(t0) > budget_s, break; end
    try
        evalc('info = os_types_hsimple_claude_auto(dsl{k});');
        fprintf(lg, '%s\tOK\tanchor %.1f\th %d\tcx %d\torth %d\tother %d\tamb %d\trf %d\n', dsl{k}, info.anchor_axis, ...
            info.n_h_simple, info.n_complex_vOS, info.n_orth_simple, info.n_other, info.ambiguous_anchor, info.has_rf);
        done = done + 1;
    catch ME
        fprintf(lg, '%s\tERROR\t%s\n', dsl{k}, ME.message);
        if ~exist([oa 'claude_types/' dsl{k}], 'dir'), mkdir([oa 'claude_types/' dsl{k}]); end
        fid = fopen([oa 'claude_types/' dsl{k} '/v2_ERROR.txt'], 'w'); fprintf(fid, '%s\n', ME.message); fclose(fid);
        done = done + 1;
    end
end
fclose(lg); fprintf('%d of %d datasets done\n', done, numel(dsl));
end
