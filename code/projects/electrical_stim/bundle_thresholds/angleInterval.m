function [angleVector] = angleInterval(inputVector)
% gives an angle vector with the same size as the input vector (degrees)
% A value of 1 in angleVector means the angle is in the 1st 2 quadrants
% with the relative axis
% A value of -1 in angleVector means the angle is in the last 2 quadrants
% with the relative axis
% A value of 0 in angleVector means the angle is 0 on the relative axis or
% NaN

    angleVector = zeros(size(inputVector));
    for i = 1:length(inputVector)
       if (isnan(inputVector(i))) 
           angleVector(i) = 0;
       elseif (inputVector(i) >= 0) && (inputVector(i) < 180) 
           angleVector(i) = 1;
       else
           angleVector(i) = -1;
       end
    end
end