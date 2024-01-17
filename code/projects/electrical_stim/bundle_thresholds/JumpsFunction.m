% Jumps function for detecting propagation
% Inputs: Jump sequence - a list of how many total jumps at each ampl
%         start - first index where TimeWindow test passed
%         maxJumps - maximum jump for this particular band

function out = JumpsFunction(JumpSequence, start, maxJumps)
if maxJumps >= 13
    candidates = find(JumpSequence < 6); % Find ampls that dont pass test
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 12
    candidates = find(JumpSequence < 6); % Find ampls that dont pass test
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 11
    candidates = find(JumpSequence < 5);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 10
    candidates = find(JumpSequence < 5); % Find ampls that dont pass test
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 9
    candidates = find(JumpSequence < 4);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 8
    candidates = find(JumpSequence < 4);
   if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 7
    candidates = find(JumpSequence < 3);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 6
    candidates = find(JumpSequence < 3);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 5
    candidates = find(JumpSequence < 2);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
elseif maxJumps >= 4
    candidates = find(JumpSequence < 2);
    if max(max(candidates)+1,start) >= length(JumpSequence)
        out = inf;
    else
        out = max(max(candidates)+1,start);  % First ampl which passes test
    end
else
    out = NaN;
end
end