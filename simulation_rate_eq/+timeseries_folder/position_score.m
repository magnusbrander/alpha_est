function [score,indSmallPosDiff] = position_score(prevPos,curPos,lenEst,originalLength)



% Assume that the orginal length correspond to a "typical" diffusion for an
% intact lambda DNA 
D_lambda = 0.1; %micrometer^2/s
% Assume the diffusion constant scales as ~1/L 
D_current = D_lambda *(originalLength/lenEst);

% Specify frame rate
frameRate = 10;
t = 1/frameRate;
% Set pixel size in same units as the diffusion constant
pixelSize = 0.16;

sigmaDif = sqrt(2*D_current*t/(pixelSize^2));

% set detection uncertainty in pixels
epsilon = 2;


% Set tollerance in movement based on length 
tollerance = ceil(3*sigmaDif + epsilon)

% Calculate difference in position between frames
posDiff = abs(prevPos - curPos)

% Fragments within tollerance
indSmallPosDiff = find(posDiff<=tollerance);

% Set logical score
if ~isempty(indSmallPosDiff)
    score = 1
else
    score = 0
end



end