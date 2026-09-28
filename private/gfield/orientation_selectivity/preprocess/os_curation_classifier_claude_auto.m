%% os_curation_classifier_claude_auto.m
% Learns GDF's earlier curation decisions (definite / maybe / rejected) from grating
% tuning features and applies them to the new datasets. Run after
% os_build_feature_table_claude_auto. Evaluation is leave-one-DATASET-out, so no
% model is ever scored on cells from a recording it was trained on.
% Two binomial logistic models on standardized features:
%   P_def: definite vs (maybe + rejected), among threshold-passing cells
%   P_os : (definite + maybe) vs rejected
% plus a 2-feature model (osi_best, corr_best: the pipeline's own metrics) and
% PCA / k-means on the aligned tuning curves as an unsupervised check.
% 2026-09-27 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/';
load([oa 'claude_features/feature_table.mat'], 'T');
O = T(~T.is_new, :);  N = T(T.is_new, :);
feats = {'osi_best','corr_best','F2n_best','F1n_best','log_rate_best','split_half_best', ...
         'dsi_max','F2n_mean','frac_cond_os','ori_consistency'};
% --- PCA on aligned, normalized tuning curves (all cells, old + new)
[coeff, score, ~, ~, expl] = pca(T.curve);
npc = 3;
T.pc = score(:, 1:npc); O.pc = T.pc(~T.is_new, :); N.pc = T.pc(T.is_new, :);
Xall = @(A) [A{:, feats}, A.pc];
fnames = [feats, arrayfun(@(k) sprintf('curvePC%d', k), 1:npc, 'UniformOutput', false)];
XO = Xall(O); XN = Xall(N);
mu = mean(XO, 1, 'omitnan'); sd = std(XO, 0, 1, 'omitnan');
z = @(X) (fillmissing(X, 'constant', 0) - mu) ./ sd;
ZO = z(XO); ZN = z(XN);
ydef = double(O.label == 2); yos = double(O.label >= 1);
ds = categorical(O.dataset); uds = categories(ds);
% --- leave-one-dataset-out predictions
Pdef = nan(height(O), 1); Pos = Pdef; P2 = Pdef;
for k = 1:numel(uds)
    te = ds == uds{k}; tr = ~te;
    b = glmfit(ZO(tr, :), ydef(tr), 'binomial'); Pdef(te) = glmval(b, ZO(te, :), 'logit');
    b = glmfit(ZO(tr, :), yos(tr), 'binomial');  Pos(te)  = glmval(b, ZO(te, :), 'logit');
    b = glmfit(ZO(tr, 1:2), ydef(tr), 'binomial'); P2(te) = glmval(b, ZO(te, 1:2), 'logit');
end
[~,~,~,R.auc_def_full] = perfcurve(ydef, Pdef, 1);
[~,~,~,R.auc_def_2feat] = perfcurve(ydef, P2, 1);
[~,~,~,R.auc_os_full] = perfcurve(yos, Pos, 1);
% conservative cutoff: smallest t with LODO precision >= 0.90 for "definite"
ts = 0.05:0.01:0.99; prec = nan(size(ts)); rec = prec;
for i = 1:numel(ts)
    c = Pdef >= ts(i); if sum(c) < 20, break; end
    prec(i) = mean(ydef(c) == 1); rec(i) = sum(ydef(c) == 1) / sum(ydef == 1);
end
it = find(prec >= 0.90, 1); R.t_def = ts(it); R.prec_def = prec(it); R.rec_def = rec(it);
% rejection cutoff: P_os below which LODO precision for "rejected" >= 0.90
prr = nan(size(ts)); rrr = prr;
for i = 1:numel(ts)
    c = Pos < ts(i); if sum(c) < 20, continue; end
    prr(i) = mean(yos(c) == 0); rrr(i) = sum(yos(c) == 0) / sum(yos == 0);
end
ir = find(prr >= 0.90, 1, 'last'); R.t_rej = ts(ir); R.prec_rej = prr(ir); R.rec_rej = rrr(ir);
% --- final models on all earlier data, applied to new candidates
bdef = glmfit(ZO, ydef, 'binomial'); bos = glmfit(ZO, yos, 'binomial');
N.P_def = glmval(bdef, ZN, 'logit'); N.P_os = glmval(bos, ZN, 'logit');
N.proposed = repmat("maybe", height(N), 1);
N.proposed(N.P_def >= R.t_def) = "definite";
N.proposed(N.P_os < R.t_rej & N.proposed ~= "definite") = "likely_not_OS";
O.P_def_lodo = Pdef; O.P_os_lodo = Pos;
R.coef_def = array2table(bdef', 'VariableNames', [{'intercept'}, fnames]);
R.coef_os  = array2table(bos',  'VariableNames', [{'intercept'}, fnames]);
R.pca_explained = expl(1:5)';
% --- unsupervised check: k-means on curve PCs, label composition per cluster
rng(1); kc = kmeans(T.pc(:, 1:npc), 4, 'Replicates', 20);
T.cluster = kc; O.cluster = kc(~T.is_new);
R.cluster_table = crosstab(O.cluster, O.label);
save([oa 'claude_features/curation_classifier_results.mat'], 'R', 'O', 'N', 'T', 'coeff', 'fnames', 'mu', 'sd', 'bdef', 'bos');
