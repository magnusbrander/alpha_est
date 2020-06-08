function bigBoxesReduced = exclude_box_overlap(firstImage,bigBoxes)

% This function deletes all overlapping boxes as well as boxes that stretch
% outside the original image.

% Check big boxes for overlap as well as out of bound in first frame
[rOrginal,cOriginal] = size(firstImage);
originalImage = [1,1,cOriginal-1,rOrginal-1]; % four-vector for the image itself


% Compute overlap between the boxes themselves. 
overlapRegions = rectint(bigBoxes,bigBoxes);

% Identify all non zero elements 
[rOverlap,cOverlap] = find(overlapRegions);

% For non diagonal terms in the outputRegions matrix we have an overlap 
indexOverlap = rOverlap(rOverlap~=cOverlap);

% Make sure we do not delete a box twice because it overlaps with several
% boxes.
indexOverlap = unique(indexOverlap);

% Delete overlapping boxes
bigBoxes(indexOverlap,:) = [];



[lengthOfBigBoxes,~] = size(bigBoxes);

fprintf('Number of big boxes after removing overlap %i.\n',lengthOfBigBoxes)


% Compute overlap region of the big boxes with the image and compare this
% to the size of the boxes. If the difference is not zero, the box is partly
% outside the image, thus remove.  
outOfImage = bigBoxes(:,3).*bigBoxes(:,4) - rectint(bigBoxes,originalImage);

% Use find to identify non zero entries
indexOutOfImage = find(outOfImage);

% Delete all boxes that stretch outside the image
bigBoxes(indexOutOfImage,:) = [];
[lengthOfBigBoxes,~] = size(bigBoxes);

fprintf('Number of big boxes after removing out of bound %i.\n',lengthOfBigBoxes)

% Return the completely non overlapping set of big boxes
bigBoxesReduced = bigBoxes; 

end


