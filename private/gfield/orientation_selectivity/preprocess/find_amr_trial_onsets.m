function [onset_idx, report] = find_amr_trial_onsets(triggers, stimulus_path, varargin)
% find_amr_trial_onsets  Robust trial-onset detection for newer-rig (540-frame)
% drifting-grating runs, with per-trial validation against the stimulus file.
%
% usage: [onset_idx, report] = find_amr_trial_onsets(triggers, stimulus_path)
%
% Background: on the newer rig the trigger line pulses once per grating
% temporal cycle (every TP/60 s at 60 Hz), so within a trial the
% inter-trigger interval is 1 s (TP 60) or ~4 s (TP 240). The older 2-s-gap
% rule (round(diff)==2) only catches boundaries whose gap happens to round
% to 2 s and silently merges trials otherwise (e.g. a truncated first trial).
% Here a trial boundary is ANY interval that is not a within-trial cycle
% interval, and each detected trial is checked against the stimulus file:
% it must contain ceil(frames/TP) triggers for the TP listed for that trial.
%
% Dropped-trigger repair: a single trigger missing INSIDE a TP-60 trial
% leaves a ~2 s within-trial gap that looks exactly like a trial boundary,
% splitting one trial into two. A detected onset j is treated as spurious
% and removed when the intervals on either side of it are both shorter than
% a normal trial and together sum to one normal trial duration (within tol).
% A genuinely truncated first trial (recording started late) does not
% satisfy this and is left alone.
%
% outputs:  onset_idx - indices into triggers of each detected trial onset
%           report    - struct: n_onsets, n_file_trials, count_mismatch_trials
%                       (empty if every trial matches; NaN if n differs),
%                       observed_counts, expected_counts, onset_intervals
%
% optional: 'tol' (default 0.1 s) tolerance for matching a cycle interval
%
% 2026-09-27 GDF + Claude

p = inputParser;
p.addParameter('tol', 0.1, @isnumeric);
p.parse(varargin{:});

triggers = triggers(:)';
txt = fileread(stimulus_path);
fr = regexp(txt, ':frames\s+(\d+)', 'tokens', 'once');
frames = str2double(fr{1});
tok = regexp(txt, ':TEMPORAL_PERIOD\s+(\d+)', 'tokens');
tp_seq = cellfun(@(c) str2double(c{1}), tok);

cycle_s = unique(tp_seq) / 60;
d = diff(triggers);
within = false(size(d));
for c = cycle_s
    within = within | abs(d - c) < p.Results.tol;
end
onset_idx = [1, find(~within) + 1];

% remove spurious onsets created by a dropped within-trial trigger
merged_at = [];
changed = true;
while changed
    changed = false;
    iv = diff(triggers(onset_idx));
    T = median(iv);
    for j = 2:numel(onset_idx)-1
        a = iv(j-1); b = iv(j);
        if a < T - p.Results.tol && b < T - p.Results.tol && abs(a + b - T) < p.Results.tol
            onset_idx(j) = [];
            merged_at(end+1) = j - 1; %#ok<AGROW>
            changed = true;
            break
        end
    end
end

report.n_onsets        = numel(onset_idx);
report.n_file_trials   = numel(tp_seq);
report.observed_counts = [diff(onset_idx), numel(triggers) - onset_idx(end) + 1];
report.expected_counts = ceil(frames ./ tp_seq);
report.onset_intervals = diff(triggers(onset_idx));
report.merged_trials   = merged_at;   % trials that absorbed a spurious split
iv = report.onset_intervals;
report.irregular_intervals = find(abs(iv - median(iv)) > p.Results.tol);
if report.n_onsets == report.n_file_trials
    report.count_mismatch_trials = find(report.observed_counts ~= report.expected_counts);
else
    report.count_mismatch_trials = NaN;
end
