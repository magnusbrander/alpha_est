function [allTimesForCutVis,allTimesForCutTrue] = extract_timeseries_from_simulation(combinedFragments)

% ######################################################################
% FUNCTION:
% Description:
% Input:
% Output:
% ######################################################################



% get total number of simulation instances
totNrInstances = 0;
for i=1:length(combinedFragments)
    % count number of saved time instances    
    totNrInstances = totNrInstances + length(combinedFragments{i}{3}); 
end


% Create vector for all times, number of visual and true number of fragments 
allTimes = zeros(1,totNrInstances);
allFragCountsVis = zeros(1,totNrInstances);
allFragCountsTrue = zeros(1,totNrInstances);

count = 1;
for i=1:length(combinedFragments)
    
    for t=1:length(combinedFragments{i}{3})
        
        allTimes(count) = combinedFragments{i}{3}(t);
        allFragCountsVis(count) = length(combinedFragments{i}{1}{t});
        allFragCountsTrue(count) = i;
        count = count + 1;
        
    end
    
    
end


allFragCountsVisCumMax = cummax(allFragCountsVis);

diffForCutVis = [0,diff(allFragCountsVisCumMax)]; 
diffForCutTrue = [0,diff(allFragCountsTrue)];

allTimesForCutVis = allTimes(diffForCutVis>0);
allTimesForCutTrue = allTimes(diffForCutTrue>0);



end