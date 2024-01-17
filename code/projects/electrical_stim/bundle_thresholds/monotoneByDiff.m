function [differenceVector] = monotoneByDiff(inputVector)
% gives a difference vector with the same size as the input vector
% A value of 1 in differenceVector means it is monotnically non-decreasing
% A value of 0 in differenceVector means it is decreasing

    differenceVector = 10*ones(size(inputVector));
    
    D = diff(inputVector);
    if (inputVector(1) > 0)
        differenceVector(1) = 1;
    end
    for i = 1:length(D)
       if (D(i) >= 0) && (inputVector(i+1) >0)
           differenceVector(i+1) = 1;
       end
    end
end