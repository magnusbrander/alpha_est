
function simulatedTimeSeries = test_script(constVal)
% Simulate cutting time distributions


% -------Constants -----------

% DNA length
L = 48500; %[basepairs]
pixelSize = 160; %[nanometers]
frayingDist = 1; %[basepairs]
diffusionConst = 1e5; %[nanometers^2 / s]
numberOfFragments = 7;
numberOfCuts = numberOfFragments - 1;

% buffer concentration
salt = 2;
k = 30;
alpha = k/(2*L);



% Generate a number of time series

nrOfTimeSeries = 100;


import timeseries_folder.*;
% myCluster = parpool('AttachedFiles','./+timeseries_folder');
% listAutoAttachedFiles(myCluster)

tic;
simulatedTimeSeries = cell(1,nrOfTimeSeries);
parfor i = 1:nrOfTimeSeries
    % simulate a single time serie
    i
    [simulatedTimeSeries{i}{1},simulatedTimeSeries{i}{2}] = generate_simulated_timeseries(alpha,...
        frayingDist,L,numberOfFragments,pixelSize,diffusionConst,salt);
    
    
end
timer = toc



end
