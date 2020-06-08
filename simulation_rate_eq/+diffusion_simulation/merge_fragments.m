function combinedFragments = merge_fragments(fragmentEvolution,resolvLim,fragLengths)
% ########################################################################
% FUNCTION: combinedFragments = merge_fragments(fragmentEvolution,resolvLim,fragLengths)
%
% Description: This function combines fragments that appears closer than
% the resolution limit to another fragment 
%
% Input: Input the fragmentEvolution cell array from diffusion_process.m,
% the resolution limit and the fragment lengths 
%
% Output: The function outputs the combined fragments and the corresponding
% time for each instance.
% ########################################################################



% Loop through all instances with multiple fragments

combinedFragments = cell(1,length(fragmentEvolution));

for p=1:length(fragmentEvolution{1}{2})
    
    % Insert positions for first fragment
    combinedFragments{1}{1}{p} = fragmentEvolution{1}{1}(p);
    % Insert lengths for first fragment
    combinedFragments{1}{2}{p} = fragLengths{1};
    % Insert times for the first fragment
    combinedFragments{1}{3} = fragmentEvolution{1}{2};
    
end


% Allocate memory in combinedFragments cell array
for i=2:length(fragmentEvolution)
    
    combinedFragments{i}{1} = cell(1,length(fragmentEvolution{i}{2}));
    combinedFragments{i}{2} = cell(1,length(fragmentEvolution{i}{2}));
    
end


import timeseries_folder.centOfMass;
for i=2:length(fragmentEvolution)
    
    currentTime = fragmentEvolution{i}{2};
    currentPositions = fragmentEvolution{i}{1};
    currentLengths = fragLengths{i};
    
    for t=1:length(currentTime)
        
        meltedPositions = currentPositions(t,:);
        meltedLengths = currentLengths;
                
        k = 1;
        while k < length(meltedPositions)
            
            
            % Compute right edge of current fragment and left edge of next
            % fragment 
            
           
            
            rightEdgeCurr = meltedPositions(k) + 0.5*meltedLengths(k);
            leftEdgeNext = meltedPositions(k+1) - 0.5*meltedLengths(k+1);
            
            if abs(rightEdgeCurr-leftEdgeNext) < resolvLim
                
                leftEdgeCurr = meltedPositions(k) - 0.5*meltedLengths(k);
                rightEdgeNext = meltedPositions(k+1) + 0.5*meltedLengths(k+1);
                newMergedLength = abs(rightEdgeNext-leftEdgeCurr);
                newMergedPos = leftEdgeCurr + 0.5*newMergedLength;
                tempPositions = zeros(1,length(meltedPositions)-1);
                tempLengths = tempPositions;
                
                count = 1;
                for j = 1:length(tempPositions)
                    
                    if j == k
                        
                        tempPositions(j) = newMergedPos; 
                        tempLengths(j) = newMergedLength; 
                        count = count + 2;
                        
                    else
                        tempPositions(j) = meltedPositions(count); 
                        tempLengths(j) = meltedLengths(count);
                        count = count + 1;
                        
                    end
                     
                end
                
                meltedPositions = tempPositions;
                meltedLengths = tempLengths;
                
                % restar the merging procedure since we might want to merge
                % the current fragment with more fragment
                k = 0;
                
                
            end
            
            k = k+1;
            
        end
        
        % Exclude fragments smaller than the resolution limit
        meltedPositions(meltedLengths < resolvLim) = [];
        meltedLengths(meltedLengths < resolvLim) = [];
        % The structure is: {number of fragment instance}{position=1/length=2/time=3}{current fragment position/lengths}
        combinedFragments{i}{1}{t} = meltedPositions;
        combinedFragments{i}{2}{t} = meltedLengths;
        combinedFragments{i}{3} = currentTime;
    end
     
    
    
end



end