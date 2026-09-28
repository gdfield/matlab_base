function T = os_build_feature_table_claude_auto(fdir, oa, new_names)
% os_build_feature_table_claude_auto  One row per threshold-passing cell, all datasets.
%
% Candidate definition = the curation pipeline's rule on find_OS_RGCs outputs:
%   any OSI > 0.30 & all DSI < 0.25 & any corr > 0.30 (over SP x TP conditions).
% Labels (earlier datasets only): 2 = definite (os_cell_list), 1 = maybe
% (os_maybe_list), 0 = rejected (passed thresholds, in neither list).
% New datasets get label NaN. Listed cells that fail the current thresholds are
% outside the candidate domain and are excluded (counted in T.Properties.UserData).
%
% Per-cell features ("best" = condition maximizing F2n * max(corr,0)):
%   osi_best, corr_best, F2n_best, F1n_best, rate_best (log10 Hz), split_half_best,
%   dsi_max, F2n_mean, frac_cond_os (fraction of conditions with OSI > 0.3),
%   ori_consistency (F2n-weighted resultant length of preferred orientation across
%   conditions; NaN when only one condition), pref_ori_best (deg, stimulus coords),
%   curve (1x12 tuning curve at best condition, normalized to max, rotated so the
%   preferred direction is first, resampled to 30-deg steps by trigonometric
%   interpolation so 8- and 12-direction data share one space).
% 2026-09-27 GDF + Claude
d = dir([fdir '*_features.mat']);
rows = {}; excl = 0;
for k = 1:numel(d)
    w = whos('-file', [fdir d(k).name]); if ~any(strcmp({w.name}, 'F')), continue; end
    S = load([fdir d(k).name], 'F'); F = S.F; ds = F.dataset;
    ncell = numel(F.cell_ids); nsp = numel(F.sps); ntp = numel(F.tps);
    osi = reshape(F.osi, ncell, []); dsi = reshape(F.dsi, ncell, []); cr = reshape(F.corr, ncell, []);
    f2 = reshape(F.F2n, ncell, []); f1 = reshape(F.F1n, ncell, []); rt = reshape(F.rate, ncell, []);
    sh = reshape(F.split_half, ncell, []); po = reshape(F.pref_ori, ncell, []);
    pass = any(osi > 0.3, 2) & all(dsi < 0.25, 2) & any(cr > 0.3, 2);
    isnew = ismember(ds, new_names);
    if ~isnew
        M = load([oa ds '_os.mat']); cf = M.os_cell_list(:)'; mb = M.os_maybe_list(:)';
        excl = excl + sum(ismember([cf mb], F.cell_ids(~pass)));
    end
    tun = reshape(F.tuning, 1, []);
    nd = numel(F.dirs); th = deg2rad(F.dirs(:)');
    for c = find(pass)'
        [~, b] = max(f2(c, :) .* max(cr(c, :), 0));
        y = double(tun{b}(c, :));
        % trigonometric interpolation to 12 points, aligned to preferred direction
        Y = fft(y); nh = floor(nd / 2);
        [~, im] = max(y); phi0 = th(im);
        g = deg2rad(0:30:330) + phi0; cv = real(Y(1)) / nd * ones(1, 12);
        for h = 1:nh
            wgt = 2; if h == nd / 2, wgt = 1; end
            cv = cv + wgt / nd * real(Y(h + 1) * exp(1i * h * (g - th(1))));
        end
        cv = cv / max(cv);
        if ncell > 0 && nsp * ntp > 1
            oc = abs(sum(f2(c, :) .* exp(2i * deg2rad(po(c, :))))) / sum(f2(c, :));
        else
            oc = NaN;
        end
        lab = NaN;
        if ~isnew
            id = F.cell_ids(c);
            if ismember(id, cf), lab = 2; elseif ismember(id, mb), lab = 1; else, lab = 0; end
        end
        rows(end+1, :) = {ds, isnew, F.cell_ids(c), lab, nd, nsp * ntp, F.reps, osi(c, b), cr(c, b), f2(c, b), f1(c, b), ...
            log10(max(rt(c, b), 0.01)), sh(c, b), max(dsi(c, :)), mean(f2(c, :)), mean(osi(c, :) > 0.3), oc, po(c, b), cv}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'dataset', 'is_new', 'cell_id', 'label', 'n_dirs', 'n_cond', 'reps', ...
    'osi_best', 'corr_best', 'F2n_best', 'F1n_best', 'log_rate_best', 'split_half_best', 'dsi_max', 'F2n_mean', ...
    'frac_cond_os', 'ori_consistency', 'pref_ori_best', 'curve'});
T.Properties.UserData.listed_cells_failing_thresholds = excl;
end
