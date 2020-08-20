function simulatedData = empirical_rate_equation(...
    alpha,constVal,statsNr)

% ########################################################################
%
% FUNCTION: simulatedData = empirical_rate_equation(alpha,constVal,statsNr)
%
% Description: 
%
% Input: 
%
% Output:  
%
% ########################################################################


% Import needed functions
import timeseries_folder.generate_sim_timeseries_lim_time_res;
import diffusion_simulation.compute_time_instances;



% Compute the discretized time interval
disTimeInt = compute_time_instances(constVal.timeResolution,constVal.maxTime);

% Arrays to store the total number of cumulative cut counts in all time intervals 
cutCountsObs = zeros(1,length(disTimeInt));
nrVisibleCuts = zeros(1,length(disTimeInt));
cutCountsTrue = zeros(1,length(disTimeInt));

% Generate many cumulative max time series
parfor n=1:statsNr
    
    [timeSeriesObs,timeSeriesTrue,visibleCuts] = generate_sim_timeseries_lim_time_res(...
        alpha,constVal);
    
    cutCountsObs(n,:) = timeSeriesObs;
    nrVisibleCuts(n,:) = visibleCuts;
    cutCountsTrue(n,:) = timeSeriesTrue;
             
end

% Compute the average cumulative max count for all times for both the observed
% and true average cut counts arrays
avgCutCountsObs = mean(cutCountsObs,1);
varCutCountsObs = var(cutCountsObs,1);
avgNrVisibleCuts = mean(nrVisibleCuts,1);
varNrVisibleCuts = var(nrVisibleCuts,1);
avgCutCountsTrue = mean(cutCountsTrue,1);
varCutCountsTrue = var(cutCountsTrue,1);



rateEqObs = gradient(avgCutCountsObs,constVal.timeResolution);
rateEqTrue = gradient(avgCutCountsTrue,constVal.timeResolution);

smoothRateEqObs = conv(rateEqObs,0.1*ones(1,10),'same');
smoothRateEqTrue = conv(rateEqTrue,0.1*ones(1,10),'same');

% Save all the simulated properties in a struct before returning
simulatedData.avgCutCountsObs = avgCutCountsObs;
simulatedData.varCutCountsObs = varCutCountsObs;
simulatedData.avgNrVisibleCuts = avgNrVisibleCuts;
simulatedData.varNrVisibleCuts = varNrVisibleCuts;
simulatedData.avgCutCountsTrue = avgCutCountsTrue;
simulatedData.varCutCountsTrue = varCutCountsTrue;
simulatedData.rateEquation = rateEqObs;
simulatedData.smoothRateEquation = smoothRateEqObs;
simulatedData.smoothRateEquationTrue = smoothRateEqTrue;
simulatedData.timeInstances = disTimeInt;

% Analytical rate
rateEqExact = 2*alpha*constVal.L*(1-exp(-alpha*constVal.frayingDist.*disTimeInt));






end