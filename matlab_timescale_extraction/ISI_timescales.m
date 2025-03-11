function [out] = ISI_timescales(timestamps)
% ISI_TIMESCALES Computes inter-spike interval (ISI) timescales and fits an exponential decay model
% Based on the algorithm introduced in Fontanier et al., 2022 (eLife)
%
% INPUT:
%   timestamps - vector of spike times (in seconds)
%
% OUTPUT:
%   out - structure containing ISI histogram, smoothed data, fit parameters, and time constant (Tau)

%- Parameters for the Exponential Fit
nPerm = 100;          % Number of random initializations for the fit
A_test = 0:0.01:20;   % Range for parameter A
B_test = 0:0.01:20;   % Range for parameter B
Tau_test = 0:0.1:1000;% Range for time constant (Tau)
maxFailure = 100;     % Maximum number of failed fits before stopping

%- Ensure timestamps are in row format
if ~isempty(timestamps)
     if size(timestamps,2) == 1
          timestamps = timestamps'; % Convert to row vector if needed
     end
end

%- Define Binning Parameters for ISI Histogram
bins_hist = [0:0.0033:1]; % Bins from 0 to 1 sec, with 3.33 ms step
FAILED = false;            % Track if fitting fails

if ~isempty(timestamps)
     %- Compute Time Differences for Lags Up to 100 Spikes
     data_all = [];
     for i = 1:100
          % Compute pairwise ISI for lags up to 100
          data_all = [data_all timestamps(i+1:end) - timestamps(1:end-i)];
     end

     % Remove ISIs that are negative or above 1 sec
     data_all = data_all(data_all > 0 & data_all < 1);

     %- Compute Histogram of ISIs
     [data, bins] = hist(data_all, bins_hist); % Histogram with defined bins

     % Min-max normalize the histogram data
     data_norm = (data - min(data)) / (max(data) - min(data));

     % Smooth the normalized data using LOESS (local regression)
     data_norm_smoothed = smooth(bins, data_norm, 0.1, 'loess'); % 10% smoothing

     %- Find Peak and Remove Initial Data Points Before Peak
     start = find(data_norm_smoothed == max(data_norm_smoothed), 1, 'first'); % Find first peak
     data_short_norm = data_norm_smoothed(start:end); % Keep data after peak
     bins_short = bins(start:end) * 1000; % Convert bins to milliseconds

     %- Define the Exponential Decay Fit Function
     % Formula: R(ΔK) = A * [exp(-Δtime / Tau) + B]
     g = fittype('A*(exp(-x/time_constant) + B)');

     %- Perform Fit with Multiple Random Initializations
     clear model_error models
     h = 0; nbFailed = 0;
     disp('    Fitting in progress')

     while h < nPerm && nbFailed < maxFailure
          try
               % Fit the model with random starting values from predefined ranges
               [f0, gof] = fit(bins_short', data_short_norm, g, ...
                    'StartPoint', [A_test(randperm(length(A_test), 1)), ...
                    Tau_test(randperm(length(Tau_test), 1)), ...
                    B_test(randperm(length(B_test), 1))]);

               h = h + 1; % Successful fit count
               model_error(h) = gof.sse; % Store sum of squared errors
               models(h).f0 = f0; % Store model parameters
               model_r(h) = gof.adjrsquare; % Store R-squared value
          catch
               nbFailed = nbFailed + 1; % Count failed fits
          end
     end

     %- Select Best Fit Model (Lowest Error)
     if nbFailed ~= maxFailure
          disp('    Fitting done')
          best_model = find(model_error == min(model_error), 1, 'first');
          f0 = models(best_model).f0;
          sse = model_error(best_model);
          r2 = model_r(best_model);

          % Display the best fit parameters
          disp(f0)

          %- keep results in output structure
          out.bins_short = bins_short;
          out.bins = bins;
          out.data_norm = data_norm;
          out.data_norm_smoothed = data_norm_smoothed;
          out.data_short_norm = data_short_norm;
          out.tau = f0.time_constant; % Extract fitted Tau
          out.lat = bins_short(1);    % First ISI peak
          out.f0 = f0;                % Best fit model
          out.sse = sse;              % Sum of squared errors
          out.r2 = r2;                % R-squared value

     else
          disp('    Fitting failed')
          FAILED = true;
     end

else
     % If timestamps are empty, set FAILED flag
     FAILED = true;
end

%- if not working, empty matrices
if FAILED
     out.bins_short = [];
     out.bins = [];
     out.data_norm = [];
     out.data_norm_smoothed = [];
     out.data_short_norm = [];
     out.tau = [];
     out.lat = [];
     out.f0 = [];
     out.sse = [];
     out.r2 = [];
end

end






