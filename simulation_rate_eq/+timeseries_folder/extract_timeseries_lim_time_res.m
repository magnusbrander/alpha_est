function [timeSeriesObsCumMax,timeSeriesTrue,obsCutTimes,nrVisCuts] = ...
    extract_timeseries_lim_time_res(allFragPos,cutInfo,timeInstances)


% #########################################################################
%
% FUNCTION: [timeSeriesObs,timeSeriesTrue] = extract_timeseries_lim_time_res(...
%    allFragPos,cutInfo,timeInstances)
%
% Description: 
%
% Input: allFragPos = the position of all fragments in all time instances,
% cutInfo = struct from nicking simulation with cut times, timeInstances =
% all time instances we are interested in
%
% Output: timeSeriesObs = observed cumulative max time series,
% timeSeriesTrue = the true cumulative max time series
%
% #########################################################################



% Arrays for observed and true time series
nrVisCuts = zeros(1,length(timeInstances));
timeSeriesTrue = zeros(1,length(timeInstances));


% Extract the number of observed fragments in all time instance
for i=1:length(timeInstances)    
    nrVisCuts(i) = length(allFragPos{i});      
end

% Compute the cumulative max of the observed number of fragments 
 timeSeriesObsCumMax = cummax(nrVisCuts) - 1; % -1 bec we count cuts and not fragments
 nrVisCuts = nrVisCuts - 1;
 % Get true cut times
 cutTimes = cutInfo.cutTimes;
 
 % Place cuts in the closest time interval for all true cuts
 for i=1:length(cutTimes)
     
     [~,index] = min(abs(cutTimes(i)-timeInstances));
     timeSeriesTrue(index) = timeSeriesTrue(index) + 1;     
 
 end
 
 % Compute the true time series by cumulative sum of cuts
 timeSeriesTrue = cumsum(timeSeriesTrue);
 
 
 % Get cut times for the observed cuts  
 obsCutTimes = zeros(1,max(timeSeriesObsCumMax));
 count = 1;
 
 prevCount = 0;
 
 for i=1:length(timeSeriesObsCumMax)
     
     currCount = timeSeriesObsCumMax(i);
     
     if prevCount < currCount         
         for j=1:(currCount -prevCount)
             
             obsCutTimes(count) = timeInstances(i);
             count = count + 1;
             
         end         
     end
     
     prevCount = currCount;
     
 end
 
 
 
end