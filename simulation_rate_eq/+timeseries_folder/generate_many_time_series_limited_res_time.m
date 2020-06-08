function [simTimeSeriesObs,simTimeSeriesTrue] = generate_many_time_series_limited_res_time(...
    alpha,nrOfTimeSeries,constVal)


% ###################################################################
%
% FUNCTION: [simTimeSeriesObs,simTimeSeriesTrue] = generate_many_time_series_limited_res_time(...
%    alpha,nrOfTimeSeries,constVal)
%
% Description: This function generates a specified number of simulated time 
% series and returns the time series in two cell arrays containing observed
% and true time series 
%
% Input: nrOfTimeSeries = the wanted number of time series, constVal =
% the needed constant values to perform the simulation
%
% Output: simTimeSeriesObs = a cell array with all observed time series
% simTimeSeriesTrue = a cell array with all true time series
%
% ###################################################################



import timeseries_folder.generate_sim_timeseries_lim_time_res;

% cell arrays to store the simulated time series in 
simTimeSeriesObs = cell(1,nrOfTimeSeries);
simTimeSeriesTrue = cell(1,nrOfTimeSeries);

parfor i = 1:nrOfTimeSeries
    % simulate a single time serie
    [~,~,~,obsCutTimes,trueCutTimes] = generate_sim_timeseries_lim_time_res(...
    alpha,constVal)
    simTimeSeriesObs{i} = obsCutTimes;
    simTimeSeriesTrue{i} = trueCutTimes;

    
    
end










end
