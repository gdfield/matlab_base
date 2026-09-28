function done = run_pref_reliability_claude_auto(budget_s)
% run_pref_reliability_claude_auto  Resumable batch of os_pref_reliability_claude_auto
% for the cells in claude_types/axis_osi_mi_cells.csv (36 kept M9 datasets), using
% best_sp / best_tp from os_cells_master_table.mat. Per-dataset results in
% claude_types/pref_reliability/<ds>.mat (table Rl). 2026-09-28 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
od = [oa 'claude_types/pref_reliability/']; if ~exist(od, 'dir'), mkdir(od); end
C = readtable([oa 'claude_types/axis_osi_mi_cells.csv'], 'TextType', 'string'); C.dataset = string(C.dataset);
S = load([oa 'claude_features/os_cells_master_table.mat'], 'OS'); OS = S.OS;
ds = unique(C.dataset); t0 = tic; done = 0;
for k = 1:numel(ds)
    f = [od char(ds(k)) '.mat']; if exist(f, 'file'), done = done + 1; continue; end
    if toc(t0) > budget_s, break; end
    c = C(C.dataset == ds(k), :); m = strcmp(OS.dataset, ds(k)); O = OS(m, :);
    [~, lo] = ismember(c.cell_id, O.cell_id);
    evalc('Rl = os_pref_reliability_claude_auto(char(ds(k)), c.cell_id, O.best_sp(lo), O.best_tp(lo));');
    save(f, 'Rl'); done = done + 1;
end
fprintf('%d of %d datasets done\n', done, numel(ds));
end
