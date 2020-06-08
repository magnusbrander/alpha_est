function [loglike,histCellMod,waitingTimeObs] = time_prob_for_alpha_sim(...
    alpha,nrOfTimeSeries,constVal,waitingTimeExp,minCutNr,timeRes)

% ###################################################################
%
% FUNCTION: [histCounts, edges] = time_prob_for_alpha_sim(alpha,constVal)
%
% Description: This function is using the whole DNA nicking and diffusion
% utilities to generate a histogram representing probabilties to see a cut
% within a certain time interval based on simulations
%
% Input: alpha = the nicking rate value for which the probability
% histograms should be simulated, constVal = a number of numerical values
% needed for the simulation such as: Length of DNA, pixel size, diffusion
% constant, maximum number of cuts, salt concentration and fraying distance
%
% Output: LOGtimeSeriesProb = the probability for observing the given time 
%series for the specified value of alpha
%
% ###################################################################


% ----------------  Generate data ------------------------- %
% First we need to generate a statistically significant number of time
% series with enough visible cuts

switch timeRes
    case 'maxRes'
        import timeseries_folder.generate_many_time_series;
        [timeSeriesObs,~] = generate_many_time_series(alpha,nrOfTimeSeries,constVal);
    case 'limitedRes'
        import timeseries_folder.generate_many_time_series_limited_res_time;
        [timeSeriesObs,~] = generate_many_time_series_limited_res_time(...
            alpha,nrOfTimeSeries,constVal);        
end





% --------------- Return trimmed array format --------------- %

% Here we convert the time series from cell to array format excluding too
% short time series from the limited resolution set

import timeseries_folder.convert_time_series_from_cell;
timeSeriesTrimArrObs = convert_time_series_from_cell(timeSeriesObs,minCutNr);

% Compute the waiting time between cuts
waitingTimeObs = diff(timeSeriesTrimArrObs,1,2);




% -------------------------- Bin data -------------------------- %

% At this point we need to create histograms for each waiting time set
% where we normalize the number of counts in each bin such that all bin values 
% sum up to one 

% First we need to create cell array to store the binning in
histCell = cell(1,minCutNr-1);

% Bin all the waiting times into counts and edges, normalized over total
% number of counts 
for n=1:minCutNr-1
    [histCell{n}{1},histCell{n}{2}] = histcounts(waitingTimeObs(:,n),...
        'Normalization','probability','BinMethod','scott');         
end



% -------------------- Compute likelihood --------------------- %

% Here we compute the log likelihood for the given time series by
% summing up the log of each individual probability obtained from the
% waiting times
import timeseries_folder.loglikelihood_from_histogram;
[loglike,histCellMod] = loglikelihood_from_histogram(histCell,waitingTimeExp);




end
















