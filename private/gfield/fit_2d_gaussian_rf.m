function [fitParams, fittedGaussian, gof] = fit_2d_gaussian_rf(rfMap)
% FIT_2D_GAUSSIAN_RF Fits a 2D Gaussian to a receptive field map
%
% Inputs:
%   rfMap - 2D matrix of sensitivity values (receptive field map)
%
% Outputs:
%   fitParams - Structure containing fitted parameters:
%               .amplitude - peak amplitude
%               .x0 - center x position
%               .y0 - center y position
%               .sigma_x - standard deviation in x
%               .sigma_y - standard deviation in y
%               .theta - rotation angle (radians)
%               .offset - baseline offset
%   fittedGaussian - 2D matrix of fitted Gaussian values
%   gof - Goodness of fit metrics (R-squared, RMSE)

    % Get dimensions of the receptive field map
    [nRows, nCols] = size(rfMap);
    
    % Create coordinate grids
    [X, Y] = meshgrid(1:nCols, 1:nRows);
    
    % Flatten the data for fitting
    xData = X(:);
    yData = Y(:);
    zData = rfMap(:);
    
    % Remove NaN values if present
    validIdx = ~isnan(zData);
    xData = xData(validIdx);
    yData = yData(validIdx);
    zData = zData(validIdx);
    
    % Initial parameter estimates
    [maxVal, maxIdx] = max(zData);
    minVal = min(zData);
    
    % Initial guesses
    amplitude0 = maxVal - minVal;
    x0_init = xData(maxIdx);
    y0_init = yData(maxIdx);
    sigma0 = min(nRows, nCols) / 4;  % Quarter of the smaller dimension
    theta0 = 0;
    offset0 = minVal;
    
    % Initial parameter vector
    % [amplitude, x0, y0, sigma_x, sigma_y, theta, offset]
    params0 = [amplitude0, x0_init, y0_init, sigma0, sigma0, theta0, offset0];
    
    % Lower and upper bounds
    lb = [0, 1, 1, 0.5, 0.5, -pi, -inf];
    ub = [inf, nCols, nRows, nCols*2, nRows*2, pi, inf];
    
    % Fitting options
    options = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 5000);
    
    % Perform the fit
    fittedParams = lsqcurvefit(@gaussian2D, params0, [xData, yData], zData, lb, ub, options);
    
    % Extract fitted parameters
    fitParams.amplitude = fittedParams(1);
    fitParams.x0 = fittedParams(2);
    fitParams.y0 = fittedParams(3);
    fitParams.sigma_x = fittedParams(4);
    fitParams.sigma_y = fittedParams(5);
    fitParams.theta = fittedParams(6);
    fitParams.offset = fittedParams(7);
    
    % Generate fitted Gaussian on full grid
    fittedGaussian = reshape(gaussian2D(fittedParams, [X(:), Y(:)]), nRows, nCols);
    
    % Calculate goodness of fit
    residuals = zData - gaussian2D(fittedParams, [xData, yData]);
    SSresid = sum(residuals.^2);
    SStotal = sum((zData - mean(zData)).^2);
    gof.rsquared = 1 - SSresid/SStotal;
    gof.rmse = sqrt(mean(residuals.^2));
    
    % Display results
    fprintf('Fitted 2D Gaussian Parameters:\n');
    fprintf('  Amplitude: %.4f\n', fitParams.amplitude);
    fprintf('  Center: (%.2f, %.2f)\n', fitParams.x0, fitParams.y0);
    fprintf('  Sigma X: %.4f\n', fitParams.sigma_x);
    fprintf('  Sigma Y: %.4f\n', fitParams.sigma_y);
    fprintf('  Rotation: %.4f radians (%.2f degrees)\n', fitParams.theta, rad2deg(fitParams.theta));
    fprintf('  Offset: %.4f\n', fitParams.offset);
    fprintf('  R-squared: %.4f\n', gof.rsquared);
    fprintf('  RMSE: %.4f\n', gof.rmse);
end

function z = gaussian2D(params, coords)
% GAUSSIAN2D Evaluates a 2D Gaussian function
%
% params: [amplitude, x0, y0, sigma_x, sigma_y, theta, offset]
% coords: [x, y] coordinates (Nx2 matrix)

    amplitude = params(1);
    x0 = params(2);
    y0 = params(3);
    sigma_x = params(4);
    sigma_y = params(5);
    theta = params(6);
    offset = params(7);
    
    x = coords(:, 1);
    y = coords(:, 2);
    
    % Rotation transformation
    x_rot = (x - x0) * cos(theta) + (y - y0) * sin(theta);
    y_rot = -(x - x0) * sin(theta) + (y - y0) * cos(theta);
    
    % 2D Gaussian equation
    z = offset + amplitude * exp(-(x_rot.^2 / (2 * sigma_x^2) + y_rot.^2 / (2 * sigma_y^2)));
end