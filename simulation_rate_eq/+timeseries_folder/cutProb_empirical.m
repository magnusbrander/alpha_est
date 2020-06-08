function [negLogP,negSmoothLogP] = cutProb_empirical(alpha,allTimeSeries,constVal,statsNr)


% #########################################################################
%
% FUNCTION: negLogP = cutProb_empirical(alpha,allTimeSeries,constVal,statsNr)
%
% Description: 
%
% Input: 
%
% Output: 
%
% #########################################################################




% First we need to generate the empirical rate equation 
import rate_function.empirical_rate_equation;

simulatedData = empirical_rate_equation(alpha,constVal,statsNr);

% Transfer all the observed cut times to indices in the time instance
% vector with discretized times 
timeInstances = simulatedData.timeInstances;
allTimeSeriesIndices = allTimeSeries;

for i=1:size(allTimeSeries,1)
    % Get one time series times
    currentTimes = allTimeSeries(i,:);
    currentIndices = currentTimes;
    for j=1:length(currentTimes)
        [~,ind] = min(abs(timeInstances-currentTimes(j)));
        currentIndices(j) = ind;
    end
    % Transfer the indices to the array where they are stored
    allTimeSeriesIndices(i,:) = currentIndices;     
    
end


% Now we need to use Michaels formula to compute the joint probability of
% observing a number of cuts at certain time 

% Intiate the log likelihood as zero to begin with 
LogP = 0;
smoothLogP = 0;

% Extract rate equation and integral of rate eqaution
rateEq = simulatedData.rateEquation;
smoothRateEq = simulatedData.smoothRateEquation;
integralRateEq = simulatedData.avgCutCountObs;

% Loop through all time series and add all contributions to the probability
for i=1:size(allTimeSeriesIndices,1)
    
    % Get current time series data    
    currTimeSeriesIndx = allTimeSeriesIndices(i,:);
    
    for j=1:length(currTimeSeriesIndx)
        
        LogP = LogP + log(rateEq(currTimeSeriesIndx(j)));
        smoothLogP = smoothLogP + log(smoothRateEq(currTimeSeriesIndx(j)));
        
    end
    % Add the last exponential term for the current time series
    LogP = LogP - integralRateEq(currTimeSeriesIndx(end));
    smoothLogP = smoothLogP - integralRateEq(currTimeSeriesIndx(end));
    
end



% Return the negative log likelihood 
negLogP = - LogP;
negSmoothLogP = - smoothLogP;


end
