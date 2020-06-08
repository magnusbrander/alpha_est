function nrOfMolPerB = boxes_in_big_boxes_simple(B,b)




overlap = rectint(b,B);   % Compute overlap of small and large boxes in this frame
normOverlap = overlap./(b(:,3).*b(:,4)); % Calculate perecentage of enclosment

nrOfMolPerB = sum(normOverlap==1,1); % Count the total number of fully enclosed boxes




end
