function boxBoundaryOfMol = molecules_to_boxes(moleculesInFrame)

% This function converts the position of the molecules to boxes by extracting 
% the max and min pixel in both the x and y direction.



% Create cell array to store box coordinates for all molecules in all
% frames. 
boxBoundaryOfMol = cell([],1);

% Extract max and min for all molecules for both x and y in all frames
for j=1:length(moleculesInFrame)
    
    tempMolInFrame = moleculesInFrame{j}; % all molecules for one frame
    boundaryOfMol = zeros(length(tempMolInFrame),4); % array for cordinates holding the box coordinates
    % Iterate through the molecules of the current frame
    for i=1:length(tempMolInFrame)
        tempMax = max(tempMolInFrame{i}); % max in x & y for current molecule
        tempMin = min(tempMolInFrame{i}); % min in x & y for current molecule
        boundaryOfMol(i,1) = tempMin(2); % x_min
        boundaryOfMol(i,2) = tempMin(1); % y_min
        boundaryOfMol(i,3) = tempMax(2)-tempMin(2); % length of box
        boundaryOfMol(i,4) = tempMax(1)-tempMin(1); % height of box
    end
    
    % Save all coordinates of the boxes for one frame.
    boxBoundaryOfMol{j,1} = boundaryOfMol; 
    
end


end 

