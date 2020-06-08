function [foundCut,row] = check_cut(fragment,frayDist,nickRow,nickCol)

% This function check if a cut happened due to a new nicking in the specific site

halfFrayDist = 0.5*frayDist;

% Determine the opposite column
if nickCol == 1
    opCol = 2;
else
    opCol = 1;    
end



% Get the total length of the fragment
lengthOfFrag = size(fragment,1);


% Set logical value to indicate if nick was found within frying distance on
% the opposite strand
foundCut = 0; 

row =[];

% Help variable to walk with
step = 0;

% Check the positive direction first 
while (step<= halfFrayDist) && (nickRow+step <= lengthOfFrag)
    
    
    if fragment(nickRow+step,opCol) == 1
        foundCut = 1;
        row = nickRow+step;
        return
    end
    
    step = step + 1;
    
end



% If nothing was found in the positve direction check the negative one
step = 0;
while (step<= halfFrayDist) && (nickRow-step > 0)
    
    
    if fragment(nickRow-step,opCol) == 1
        foundCut = 1;
        row = nickRow-step;
        return
    end
    
    step = step + 1;
    
end



end
