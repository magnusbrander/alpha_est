function [allFragPosMerged,allFragLengthsMerged] = merge_fragments_lim_res(...
    allFragPos,allFragLengths,resolvLim)



% ########################################################################
%
% FUNCTION: [allFragPosMerged,allFragLengthsMerged] = merge_fragments_lim_res(...
%    allFragPos,allFragLengths,timeInstances)
%
% Description:
%
% Input: allFragPos = cell array with all fragment positions at all times,
% allFragLengths = cell array with all fragment lengths at all times,
% resolvLim = smallest distance that we can separate fragment on
%
% Output: allFragPosMerged = the position of the fragments after merge,
% allFragLengthsMerged = the lengths of the resulting merged fragments
%
% ########################################################################


allFragPosMerged = allFragPos;
allFragLengthsMerged = allFragLengths;

for t=1:length(allFragPos)
    
    
    
    mergedPos = allFragPos{t};
    mergedLengths = allFragLengths{t};
    
    k = 1;
    while k < length(mergedPos)
        
        % Compute right edge of current fragment and left edge of next
        % fragment
        rightEdgeCurr = mergedPos(k) + 0.5*mergedLengths(k);
        leftEdgeNext = mergedPos(k+1) - 0.5*mergedLengths(k+1);
        
        
        if abs(rightEdgeCurr-leftEdgeNext) < resolvLim
            
            leftEdgeCurr = mergedPos(k) - 0.5*mergedLengths(k);
            rightEdgeNext = mergedPos(k+1) + 0.5*mergedLengths(k+1);
            newMergedLength = abs(rightEdgeNext-leftEdgeCurr);
            newMergedPos = leftEdgeCurr + 0.5*newMergedLength;
            tempPositions = zeros(1,length(mergedPos)-1);
            tempLengths = tempPositions;
            
            count = 1;
            for j = 1:length(tempPositions)
                
                if j == k
                    
                    tempPositions(j) = newMergedPos;
                    tempLengths(j) = newMergedLength;
                    count = count + 2;
                    
                else
                    tempPositions(j) = mergedPos(count);
                    tempLengths(j) = mergedLengths(count);
                    count = count + 1;
                    
                end
                
            end
            
            mergedPos = tempPositions;
            mergedLengths = tempLengths;
            
            % restar the merging procedure since we might want to merge
            % the current fragment with more fragment
            k = 0;
            
            
        end
        
        % Move on to merging next fragment
        k = k+1;
        
    end
    
    % Exclude fragments smaller than the resolution limit
    mergedPos(mergedLengths < resolvLim) = [];
    mergedLengths(mergedLengths < resolvLim) = [];
    
    % Save the merged postions and lengths
    allFragPosMerged{t} = mergedPos;
    allFragLengthsMerged{t} = mergedLengths;
    
end




end
