function [countsInB,moleculeInBIndex] = boxes_in_bigger_boxes(b,B)


% Determine the number of big boxes
[lengthOfB,~] = size(B);
% Alocate array for the counts in big boxes for all frames
countsInB = zeros(length(b),lengthOfB);
moleculeInBIndex = cell([],lengthOfB);

for i=1:length(b)
    bTemp = b{i}; % Look at one frame at the time
    
    overlap = rectint(bTemp,B);   % Compute overlap of small and large boxes in this frame
    normOverlap = overlap./(bTemp(:,3).*bTemp(:,4)); % Normalize the output
    
    nrOfMolPerB = sum(normOverlap==1,1); % Count the total number of fully enclosed boxes
    countsInB(i,:) = nrOfMolPerB;   % Store the number of smaller boxes in the big boxes
    
    
    % Following is done in order to relate the specific molecules in each frame
    % to their right big box
    
    % Set all values for boxes that are not fully enclosed to zero
    normOverlap(normOverlap~=1) = 0;
    
    % Go through all bigger boxes to check for smaller boxes and their
    % index
    for k=1:lengthOfB
        row= find(normOverlap(:,k)); % Finds the fully enclosed smaller boxes in big box k       
        moleculeInBIndex{i,k} = row; % Save the obtained indexes
    end
end




end

