function boundaryOfMol = molecules_to_boxes_singel(moleculesInFrame)

% This function converts the position of the molecules to boxes by extracting 
% the max and min pixel in both the x and y direction.




    boundaryOfMol = zeros(length(moleculesInFrame),4); % array for cordinates holding the box coordinates
    % Iterate through the molecules of the current frame
    for i=1:length(moleculesInFrame)
        tempMax = max(moleculesInFrame{i}); % max in x & y for current molecule
        tempMin = min(moleculesInFrame{i}); % min in x & y for current molecule
        boundaryOfMol(i,1) = tempMin(2); % x_min
        boundaryOfMol(i,2) = tempMin(1); % y_min
        boundaryOfMol(i,3) = tempMax(2)-tempMin(2); % length of box
        boundaryOfMol(i,4) = tempMax(1)-tempMin(1); % height of box
    end



end 
