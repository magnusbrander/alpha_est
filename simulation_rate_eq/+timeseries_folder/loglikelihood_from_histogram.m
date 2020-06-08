function [loglikelihood,histCellMod] = loglikelihood_from_histogram(histCell,waitingTimeExp)




% ###################################################################
%
% FUNCTION: output = loglikelihood_from_histogram(input)
%
% Description: This function computes the log likelihood of the given
% experimental waiting time distrubutions from the simulated probability
% histogram 
%
% Input: nrOfTimeSeries = the wanted number of time series, constVal =
% the needed constant values to perform the simulation
%
% Output: simTimeSeriesObs = a cell array with all observed time series
% simTimeSeriesTrue = a cell array with all true time series
%
% ###################################################################


histCellMod = histCell;

loglikelihood = 0;
for i = 1:length(histCell)
    
    tempEdges = histCell{i}{2};
    binWidth = (tempEdges(2)-tempEdges(1));
    tempBinCent = 0.5*binWidth+tempEdges(1:end-1);
    tempWaitingtTimes = waitingTimeExp(:,i);
    probInBin = histCell{i}{1};
    
    
    % Rescale bin counts to exclude zero values and insert the smallest non
    % zero value instead 
    smallestNonZeroVal = min(probInBin(probInBin~=0)); % find zeros 
    probInBin(probInBin == 0) = smallestNonZeroVal;    % replace zeros
    probInBin = probInBin/(sum(probInBin)); % Normalize 
    histCellMod{i}{1} = probInBin;
     
    for n=1:length(tempWaitingtTimes)
        
        % Compute closest bin to current waiting time
        [val,index] = min(abs(tempBinCent-tempWaitingtTimes(n)));
        
        % Compute the loglikelihood and add to sum
        loglikelihood = loglikelihood + log(probInBin(index));
        
        
    end
    
   
end








end
