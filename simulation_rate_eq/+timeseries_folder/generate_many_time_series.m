
function [simTimeSeriesObs,simTimeSeriesTrue] = generate_many_time_series(alpha,nrOfTimeSeries,constVal)



% ###################################################################
%
% FUNCTION: [simTimeSeriesObs,simTimeSeriesTrue] = generate_many_time_series(...
% nrOfTimeSeries,constVal)
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



% ------- Extract constants -----------


% DNA length in basepairs 
L = constVal.L; %[basepairs]
pixelSize = constVal.pixelSize; %[nanometers]
frayingDist = constVal.frayingDist; %[basepairs]
diffusionConst = constVal.diffConst; %[nanometers^2 /s]
numberOfCuts = constVal.nrOfCuts;
salt = constVal.salt;
n0 = 2*L*constVal.intialNickDensity;


import timeseries_folder.generate_simulated_timeseries;

% cell arrays to store the simulated time series in 
simTimeSeriesObs = cell(1,nrOfTimeSeries);
simTimeSeriesTrue = cell(1,nrOfTimeSeries);

parfor i = 1:nrOfTimeSeries
    % simulate a single time serie
    [simTimeSeriesObs{i},simTimeSeriesTrue{i}] = generate_simulated_timeseries(alpha,...
        frayingDist,L,numberOfCuts,pixelSize,diffusionConst,salt,n0);
    
    
end



end
