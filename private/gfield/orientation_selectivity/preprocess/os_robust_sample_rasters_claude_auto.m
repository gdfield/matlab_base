function S = os_robust_sample_rasters_claude_auto(varargin)
% os_robust_sample_rasters_claude_auto  Raster + polar plots (plot_direction_tuning,
% best condition by r2_lbx) for a random, stratified sample of "background" cells
% (earlier datasets, never listed by GDF and failing the find_OS_RGCs gate) that the
% candidate M10 rule flags: >= 1 significant condition (q < 0.05, r2 > r1, >= 50
% spikes) AND r2_lbx >= x (0.16). Strata of r2_lbx: [0.16 0.2), [0.2 0.3), >= 0.3;
% n_per stratum (10), seed (1). Titles give r2_lbx, reliability, spikes, and the
% gate values (max OSI, max DSI, max corr over conditions) so the gate failure is
% visible. Writes claude_figures/robust/flagged_background_sample/<ds>_cell<id>.pdf,
% and sample_index.csv (page = random viewing order, stratum hidden in the order).
% 2026-09-28 GDF + Claude
p = inputParser; p.addParameter('x', 0.16); p.addParameter('n_per', 10); p.addParameter('seed', 1);
p.addParameter('edges', [0.16 0.2 0.3 Inf]); p.parse(varargin{:}); o = p.Results;
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
od = [oa 'claude_figures/robust/flagged_background_sample/']; if ~exist(od, 'dir'), mkdir(od); end
L = load([oa 'claude_robust/robust_cells.mat'], 'P'); P = L.P;
fl = find(P.label == "background" & P.k_sig >= 1 & P.r2_lbx_best >= o.x);
rng(o.seed, 'twister'); pick = [];
for s = 1:numel(o.edges) - 1
    c = fl(P.r2_lbx_best(fl) >= o.edges(s) & P.r2_lbx_best(fl) < o.edges(s + 1));
    pick = [pick; c(randperm(numel(c), min(o.n_per, numel(c))))]; %#ok<AGROW>
end
S = P(pick, {'dataset', 'cell_id', 'n_dirs', 'n_cond', 'k_sig', 'r2_best', 'r2_lbx_best', 'r1_best', 'rel_best', 'nspk_best', 'best_sp', 'best_tp', 'pref_axis_best'});
S.stratum = discretize(S.r2_lbx_best, o.edges);
S.gate_osi = nan(height(S), 1); S.gate_dsi = S.gate_osi; S.gate_corr = S.gate_osi; S.page = nan(height(S), 1);
reg = load([oa 'claude_features/registry_used.mat']);
ord = randperm(height(S)); S.page(ord) = 1:height(S);
[ud, ~, g] = unique(S.dataset);
fh = figure('Visible', 'off');
files = strings(height(S), 1);
for k = 1:numel(ud)
    dl = reg.allds(strcmp(reg.nm, ud(k)));
    evalc('dr = load_data(dl.grating_datapath); dr = load_neurons(dr); dr = os_trigger_check(dr, dl.stimulus_path, dl.trigger_interval);');
    F = load([oa 'claude_features/' char(ud(k)) '_features.mat'], 'F'); F = F.F;
    for i = find(g == k)'
        ci = find(F.cell_ids == S.cell_id(i));
        S.gate_osi(i) = max(F.osi(ci, :), [], 'all'); S.gate_dsi(i) = max(F.dsi(ci, :), [], 'all'); S.gate_corr(i) = max(F.corr(ci, :), [], 'all');
        ri = find(dr.cell_ids == S.cell_id(i));
        [tun, sn] = get_direction_tuning(dr.spikes{ri}, dr.stimulus, 'SP', S.best_sp(i), 'TP', S.best_tp(i));
        ttl = sprintf('%s cell %d  SP %g TP %g | r2lbx %.2f r2 %.2f r1 %.2f rel %.2f spikes %d sig %d/%d | gate: OSI %.2f DSI %.2f corr %.2f', ...
            S.dataset(i), S.cell_id(i), S.best_sp(i), S.best_tp(i), S.r2_lbx_best(i), S.r2_best(i), S.r1_best(i), S.rel_best(i), ...
            S.nspk_best(i), S.k_sig(i), S.n_cond(i), S.gate_osi(i), S.gate_dsi(i), S.gate_corr(i));
        plot_direction_tuning(tun, sn, dr.stimulus, 'fig_num', fh, 'fig_title', ttl, 'clear_fig', true); drawnow;
        files(i) = sprintf('%s%s_cell%d.pdf', od, S.dataset(i), S.cell_id(i));
        exportgraphics(fh, files(i), 'ContentType', 'image', 'Resolution', 150);
    end
end
close(fh);
S.file = files; writetable(S, [od 'sample_index.csv']);
end
