function newFragPos = move_frag_if_possible(prevFragPos,currLengths,...
    leftOrRight,fragToMove,jumpDist)

% ########################################################################
%
% FUNCTION: movedFragPos = move_frag_if_possible(input)
%
% Description: 
%
% Input: prevPos = , currLengths = 
%
% Output: newFragPos = position of all fragments after attempt to move
%
% ########################################################################


newFragPos = prevFragPos;


if leftOrRight == 1  % Attempt to move fragment to the left
    
    if fragToMove == 1 % If we look at left most fragment we can always move left
        
        newFragPos(fragToMove) = prevFragPos(fragToMove) - jumpDist;
        
    else % Check if it will overlap or cross other fragment
        
        % Extract current lengths and positions
        l_1 = currLengths(fragToMove-1);
        l_2 = currLengths(fragToMove);
        x1 = prevFragPos(fragToMove-1);
        x2 = prevFragPos(fragToMove);
        
        % Compute overlap and crossing checks
        
        overlapCheck = abs((x2 - 0.5*l_2) - (x1 + 0.5*l_1));  % Should be larger than jumping dist
        crossCheck = x2 - x1;  % Should be larger than jumping distance
        
        % If both checks return positive values we can move the
        % fragment to the left
        if overlapCheck > jumpDist && crossCheck > jumpDist
            newFragPos(fragToMove) = prevFragPos(fragToMove) - jumpDist;
        end
        
        
    end
    
    
else  % Attempt to move fragment to the right
    
    
    if fragToMove == length(currLengths) % If we look at right most fragment we can always move right
        
        newFragPos(fragToMove) = prevFragPos(fragToMove) + jumpDist;
        
    else % Check if it will overlap or cross other fragment
        
        % Check overlap
        l_1 = currLengths(fragToMove);
        l_2 = currLengths(fragToMove+1);
        x1 = prevFragPos(fragToMove);
        x2 = prevFragPos(fragToMove+1);
        overlapCheck = abs((x2 - 0.5*l_2) - (x1 + 0.5*l_1)); % Should be larger than jumping dist
        crossCheck = x2 - x1;  % Should be larger than zero
        
        % If both checks return positive values we can move the
        % fragment to the right
        if overlapCheck > jumpDist && crossCheck > jumpDist
            newFragPos(fragToMove) = prevFragPos(fragToMove) + jumpDist;
        end
        
    end
    
end





end