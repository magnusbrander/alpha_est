function timeInstances = compute_time_instances(dt,maxTime,startTime)

% ########################################################################
%
% FUNCTION: timeInstances = compute_time_instances(timeResolution,maxTime)
%
% Description: 
%
% Input: timeResolution = the wanted resolution in time, maxTime = the last time
%
% Output: timeInstances = the discretized time interval with resolution
% equal to timeResolution and maximum value of maxTime
%
% ########################################################################


timeInstances = (startTime+dt):dt:maxTime;


% n = length(varargin);
% switch n
%     
%     case 2
%         
%         timeResolution = varargin{1};
%         maxTime = varargin{2};
%         startTime = 0.5*timeResolution;
%         
%     case 3
%         
%         timeResolution = varargin{1};
%         maxTime = varargin{2};
%         startTime = varargin{3}+timeResolution;
%         
% end
% 
% 
% 
% t = startTime;
% count = 1;
% 
% % Iterate through time interval to know how many instances there are
% while t <= maxTime - timeResolution
%     
%     t = t + timeResolution;
%     count = count + 1;   
%     
%     
% end
% 
% % Create vector for all instances
% timeInstances = zeros(1,count);
% timeInstances(1) = startTime;
% 
% % Fill the timeInstance vector with all the time instance
% if count>1
%     for n=2:count
%         
%         timeInstances(n) = timeInstances(n-1) + timeResolution;
%                 
%     end
% end
% 
% % If the time resolution and max time does not allow for an evenly spaced
% % interval we modify the length of the last interval
% 
% if timeInstances(end) < maxTime
%     timeInstances(end+1) = maxTime;
% end
% 
% 


end