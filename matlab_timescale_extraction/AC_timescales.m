function [out] = AC_timescales(spikeCounts, binSize)

% AC_TIMESCALES calculates the temporal autocorrelation of spike counts 
% based on the method from Murray et al. 2014 (Nature Neuroscience)
% and fits an exponential decay to the autocorrelation function.
%
% INPUTS:
%   spikeCounts - A matrix of spike counts (numTrials x numTimeBins), where each row
%                 corresponds to a trial, and each column represents a time bin.
%   binSize     - The duration of each time bin in milliseconds (e.g., 50 ms).
%
% OUTPUT:
%   out - A structure containing:
%         bins (time lags in ms), data_norm (autocorrelation values), 
%         tau (decay time constant), lat (latency of decay start), 
%         f0 (fitted model), sse (sum of squared errors), r2 (fit quality).

FAILED = false; % Flag to indicate whether the process fails

%- Parameters for the Exponential Fit 
nPerm = 100;           % Number of random initializations for fitting
A_test = 0:0.01:20;    % Range for parameter A (scaling factor)
B_test = 0:0.01:20;    % Range for parameter B (offset)
Tau_test = 0:0.1:1000; % Range for time constant (Tau)
maxFailure = 100;      % Maximum number of fit failures before stopping

%- Determine Number of Time Bins 
[~, numTimeBins] = size(spikeCounts); % Extract number of time bins

% Compute the maximum lag for which autocorrelation is calculated
maxLag = numTimeBins - 1; 
R = zeros(maxLag, 1); % Preallocate autocorrelation values

%- Compute Autocorrelation Function 
for k = 1:maxLag
    % Define indices for time bins separated by lag k
    ind1 = 1:(numTimeBins - k);     % Current time bins
    ind2 = (1 + k):numTimeBins;     % Future time bins with lag k
    
    % Extract spike counts at these time bins
    N1 = spikeCounts(:, ind1); % Spike counts at initial time bins
    N2 = spikeCounts(:, ind2); % Spike counts at lagged time bins

    % Compute mean spike count across trials
    meanN1 = mean(N1, 1); % Mean for N1
    meanN2 = mean(N2, 1); % Mean for N2

    % Compute covariance and variance across trials
    covMatrix = mean((N1 - meanN1) .* (N2 - meanN2), 1);
    varN1 = mean((N1 - meanN1).^2, 1);
    varN2 = mean((N2 - meanN2).^2, 1);

    % Compute Pearson correlation coefficient (autocorrelation)
    R(k) = mean(covMatrix ./ sqrt(varN1 .* varN2));
end

%- Identify the Start of the Exponential Decay 
% Compute the difference between consecutive autocorrelation values
diffR = diff(R);
firstDecreaseIdx = find(diffR < 0, 1); % Find the first point where correlation decreases

if isempty(firstDecreaseIdx)
    firstDecreaseIdx = 1; % If no decrease found, start at the first bin
end

%- Prepare Data for Exponential Fit 
fitLags = (firstDecreaseIdx:maxLag) * binSize; % Convert lags to milliseconds
fitData = R(firstDecreaseIdx:end); % Extract autocorrelation data after first decrease

% Define the exponential decay function:
% R(Δt) = A * exp(-Δt / Tau) + B
g = fittype('A*(exp(-x/time_constant)+B)');

%-  Perform Multiple Fits and Select the Best One 
clear model_error models
h = 0; nbFailed = 0;
disp('    Fitting in progress')

while h < nPerm && nbFailed < maxFailure
    try
        % Fit the model using randomly chosen starting values
        [f0, gof] = fit(fitLags', fitData, g, ...
            'StartPoint', [A_test(randperm(length(A_test), 1)), ...
                           Tau_test(randperm(length(Tau_test), 1)), ...
                           B_test(randperm(length(B_test), 1))]); 
                       
        h = h + 1; % Increment successful fit counter
        model_error(h) = gof.sse; % Store sum of squared errors
        models(h).f0 = f0; % Store model parameters
        model_r(h) = gof.adjrsquare; % Store R-squared value (fit quality)
    catch
        nbFailed = nbFailed + 1; % Count failures
    end
end

%% --- Select Best Fit Model (Lowest Error) ---
if nbFailed ~= maxFailure
    disp('    Fitting done')

    % Find the best model based on lowest error (SSE)
    best_model = find(model_error == min(model_error), 1, 'first');
    f0 = models(best_model).f0; % Extract best-fitting model
    sse = model_error(best_model); % Extract best fit SSE
    r2 = model_r(best_model); % Extract R² value for best fit

    disp(f0) % Display best fit parameters

    %- Store Results in Output Structure 
    out.bins = binSize:binSize:maxLag * binSize; % Time lags in ms
    out.data_norm = single(R); % Autocorrelation values
    out.tau = f0.time_constant; % Extract fitted Tau (decay time constant)
    out.lat = out.bins(firstDecreaseIdx); % Latency of first drop in autocorrelation
    out.f0 = f0; % Store best fit model
    out.sse = sse; % Store SSE of best fit
    out.r2 = r2; % Store R² value of best fit

else
    disp('    Fitting failed')
    FAILED = true;
end

%- empty if failed
if FAILED
    out.bins = [];
    out.data_norm = [];
    out.tau = [];
    out.lat = [];
    out.f0 = [];
    out.sse = [];
    out.r2 = [];
end

end
