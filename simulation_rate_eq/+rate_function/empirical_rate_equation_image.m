function simulatedData = empirical_rate_equation_image(...
    alpha,constVal,statsNr)

% ########################################################################
%
% FUNCTION: simulatedData = empirical_rate_equation_image(alpha,constVal,statsNr)
%
% Description: 
%
% Input: 
%
% Output:  
%
% ########################################################################


% Import needed functions
import timeseries_folder.generate_sim_timeseries_image;
import diffusion_simulation.compute_time_instances;




% Compute the discretized time interval
startTime = 0;
disTimeInt = compute_time_instances(constVal.timeResolution,constVal.maxTime,startTime);

% Arrays to store the total number of cumulative cut counts in all time intervals 
cutCountsObs = zeros(1,length(disTimeInt));
nrVisibleCuts = zeros(1,length(disTimeInt));
cutCountsTrue = zeros(1,length(disTimeInt));

% Generate many cumulative max time series
parfor n=1:statsNr
    n
    [timeSeriesObs,timeSeriesTrue,visibleCuts] = generate_sim_timeseries_image(...
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

smoothRateEqObs = smoothdata(rateEqObs,'movmean',10);
smoothRateEqTrue = smoothdata(rateEqTrue,'movmean',10);

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
simulatedData.cutCountsTrue = cutCountsTrue;
simulatedData.cutCountsObs = cutCountsObs;
simulatedData.nrVisibleCuts = nrVisibleCuts;

% Analytical rate
rateEqExact = 2*alpha*constVal.L*(1-exp(-alpha*constVal.frayingDist.*disTimeInt));


figure; hold on; box on;
plot(disTimeInt,rateEqObs,'k.-');
plot(disTimeInt,rateEqTrue,'r.-');
plot(disTimeInt,smoothRateEqObs,'g--')
plot(disTimeInt,smoothRateEqTrue,'b--')
plot(disTimeInt,rateEqExact,'m--')
xlabel('time')
ylabel('r(t)')
legend('raw obs','raw true','smooth obs','smooth true','analytical')
title(strcat('alpha = ', num2str(alpha)));






end