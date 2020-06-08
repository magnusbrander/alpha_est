function instanceIndex = get_instance(t,timeInstances,instCount)

% ########################################################################
%
% FUNCTION: instanceIndex = get_instance(t,timeInstances)
%
% Description: 
%
% Input: t = , timeInstances = 
%
% Output: instanceIndex = 
%
% ########################################################################


instanceIndex = instCount;

for i=instCount:1:length(timeInstances)
    
    % If the time is larger or eq. than current time instance we know we are 
    % at least at the current time instance 
    if timeInstances(i) <= t
        instanceIndex = i;
    else 
        % If time is smaller than current time instance we know that we are
        % at previous instance 
        return;
    end
    
    
end


end
