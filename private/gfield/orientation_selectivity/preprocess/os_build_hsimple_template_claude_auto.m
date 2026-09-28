function tmpl = os_build_hsimple_template_claude_auto(member_ids)
% os_build_hsimple_template_claude_auto  Build the h-simple RF template from
% 2012-10-31-0 (used by os_hsimple_score_claude_auto / os_types_hsimple v3).
% member_ids: grating-run cell ids of the template cells. The template in use
% (2026-09-28) was built from the 57 v2 h-simple cells of 2012-10-31-0 that had a
% WN RF and a Vision fit (GDF-confirmed classes); their ids are stored in
% claude_types/hsimple_template_2012-10-31-0.mat (tmpl.member_ids), so calling
%   S = load(<that file>); os_build_hsimple_template_claude_auto(S.tmpl.member_ids)
% reproduces it. Each member: polarity-signed robust-z RF centered on its Vision
% fit (row = H - y + 0.5, col = x + 0.5), resampled to the members' median fit
% size, rotated so the preferred axis is at 0 deg (y up), unit norm; template =
% mean. Does not overwrite the saved template; returns it.
% 2026-09-28 GDF + Claude
oa = '/Users/gfield/Development/matlab_base/private/gfield/orientation_selectivity/os_analysis/claude_types/';
L = load([oa '2012-10-31-0/types.mat'], 'T', 'rfimg'); T = L.T; Z = L.rfimg;
h = find(ismember(T.cell_id, member_ids) & ~cellfun(@isempty, Z) & ~isnan(T.fit_x));
sz = sqrt(T.fit_sd1 .* T.fit_sd2); ref = median(sz(h));
A = nan(25, 25, numel(h));
for j = 1:numel(h)
    i = h(j); H = size(Z{i}, 1);
    A(:, :, j) = os_rf_align_claude_auto(Z{i}, T.polarity(i), H - T.fit_y(i) + 0.5, T.fit_x(i) + 0.5, T.pref_axis(i), 12, sz(i) / ref);
end
tmpl = struct('M', mean(A, 3), 'A_members', A, 'member_ids', T.cell_id(h), 'dataset', {repmat({'2012-10-31-0'}, numel(h), 1)}, ...
    'ref_size', ref, 'source', '2012-10-31-0 h_simple members; see header', 'date', char(datetime('now', 'Format', 'yyyy-MM-dd')));
end
