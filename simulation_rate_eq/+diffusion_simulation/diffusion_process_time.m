
function fragmentEvolution = diffusion_process_time(cutInfo,maxTime,intialLength,pixelSize,diffusionConst,salt)
% ########################################################################
%
% FUNCTION: fragmentEvolution = diffusion_process(cutInfo,intialLength,pixelSize,diffusionConst)  
% Input: cutInfo = class contaning cut times, cut positions cutted fragment
% number, number of nicks at cut and the complete fragment with marked
% nicks after the cut happened. initialLength = the intial length in
% basepairs, pixelSize = pixel size in nanometers, diffusionConst = empirical
% diffusion constant for a standard lambda DNA in nanometers^2/s
% Output: 
% Description:
%
% ########################################################################


% Set the diffusion jumping distance equal to the pixel size
jumpDist = pixelSize;

% Strategy: simulate intial fragment up until time for first cut, then
% simulate the diffusion for the resulting two fragments up until time for
% second cut etc... 

% extract all lengths of all fragments at each time interval and compute
% the real length in nanometers
import diffusion_simulation.extract_frag_lengths;
fragLengths = extract_frag_lengths(cutInfo.fragmentCollection,intialLength,salt);


% Intiate loop for simulating each time interval 

if cutInfo.cutTimes == 0
    cutTimes = maxTime;
else
    cutTimes = [cutInfo.cutTimes,maxTime];
end



% Create cell array for storing all positions at all times
fragmentEvolution = cell(1,length(cutTimes));

% Set time counter
t = 0;

import diffusion_simulation.assign_new_pos;
import diffusion_simulation.get_rates;
for i = 1:length(cutTimes)
    
    % calculate the rates for jumping left and right respectivly for each
    % fragment    
    [hopRatesLeft,hopRatesRight] = get_rates(fragLengths{i},diffusionConst,fragLengths{1},pixelSize);    
    mu = sum(hopRatesLeft) + sum(hopRatesRight);
    
    % set maximum simulation time up until next cut
    tMax = cutTimes(i);

    
    % Assign positions of all the fragments after cut
    
    if i == 1  % Special case if we are at intial fragment
        
        % create time and postion vector to store evolution in
        fragPos = 0;
        theTime = 0;
        currFragLengths = fragLengths{i};
        
    else % If not in first instance we need to calculate new fragment positions 
        
        % Set time for current simulation round 
        theTime = t;
        
        % Save previous last positions
        prevPositions = fragPos(end,:);
        
        % Extract previous last lengths and current lenghts
        prevLengths = fragLengths{i-1};
        currFragLengths = fragLengths{i};
        
        % Extract the number for the cutted fragment
        cuttedFragNr = cutInfo.cuttedFragmentNr(i-1);
        
        % Determine new positions
        fragPos = assign_new_pos(prevLengths,currFragLengths,prevPositions,cuttedFragNr);
        
    end

   
    while t<tMax 
        
        % Draw a time until next move using gillespie
        t = t - log(rand)/mu;
        % Save all times at which something happens
        theTime(end+1,1) = t;
        
    end
    
    if length(theTime) > 2e6
        fragmentEvolution = cell.empty;
        fprintf('Too long time array')
        return
    end
    
    
    % Allocate and compute necessary things
    % Since left and right hop rates are equal we can sample the fragment first
    % and then left or right direction as long as we get the right mu for
    % the waiting time distribution
    theFragToMoveVec = randsample(length(hopRatesLeft),length(theTime),'true',hopRatesLeft);
    leftOrRightVec = randi(2,1,length(theTime));
    fragPosAll = zeros(length(theTime)-1,length(fragPos));
    fragPos = [fragPos;fragPosAll];
    
    % Loop through all time where something happens
    for move = 2:length(theTime)
        
        
        theFragToMove = theFragToMoveVec(move);
        
        % Sample if fragment should be moved to the left or right 
        leftOrRight = leftOrRightVec(move);
        
        % Copy previous positions to the current positions
        fragPos(move,:) = fragPos(move-1,:);
        
       
        if leftOrRight == 1
            % Move fragment to the left if possible 
            
            if theFragToMove == 1 % If we look at left most fragment we can always move left
                
                fragPos(move,theFragToMove) = fragPos(move,theFragToMove) - jumpDist;
                
            else % Check if it will overlap or cross other fragment
                
                % Extract current lengths and positions
                l_1 = currFragLengths(theFragToMove-1);
                l_2 = currFragLengths(theFragToMove);
                x1 = fragPos(move,theFragToMove-1);
                x2 = fragPos(move,theFragToMove);
                
                % Compute overlap and crossing checks 
                
                overlapCheck = abs((x2 - 0.5*l_2) - (x1 + 0.5*l_1));  % Should be larger than jumping dist 
                crossCheck = x2 - x1;  % Should be larger than jumping distance 
                
                % If both checks return positive values we can move the
                % fragment to the left
                if overlapCheck > jumpDist && crossCheck > jumpDist
                    fragPos(move,theFragToMove) = fragPos(move,theFragToMove) - jumpDist;
                end
                
            end
            
   
        else
            
            % Attempt to move fragment to the right
            
            if theFragToMove == length(currFragLengths) % If we look at right most fragment we can always move right
                
                fragPos(move,theFragToMove) = fragPos(move,theFragToMove) + jumpDist;
                
            else % Check if it will overlap or cross other fragment
                
                % Check overlap
                l_1 = currFragLengths(theFragToMove);
                l_2 = currFragLengths(theFragToMove+1);
                x1 = fragPos(move,theFragToMove);
                x2 = fragPos(move,theFragToMove+1);
                overlapCheck = abs((x2 - 0.5*l_2) - (x1 + 0.5*l_1)); % Should be larger than jumping dist
                crossCheck = x2 - x1;  % Should be larger than zero 
                
                % If both checks return positive values we can move the
                % fragment to the right
                if overlapCheck > jumpDist && crossCheck > jumpDist
                    fragPos(move,theFragToMove) = fragPos(move,theFragToMove) + jumpDist;
                end
                
            end
            
        end
              
        
    end
    
    fragmentEvolution{i}{1} = fragPos; % Position difference in nanometers rela. to intial pos of 1st frag
    fragmentEvolution{i}{2} = theTime; % Time in seconds from start
    fragmentEvolution{i}{3} = [theTime,fragPos];
    
end




end

