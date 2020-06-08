function [allFragPos,allFragLengths,timeInstances] = diffusion_process_limited_res_time(cutInfo...
    ,maxTime,timeResolution,intialLength,pixelSize,diffusionConst,salt)


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

% Import needed functions 
import diffusion_simulation.extract_frag_lengths;
import diffusion_simulation.assign_new_pos;
import diffusion_simulation.get_rates;
import diffusion_simulation.move_frag_if_possible;
import diffusion_simulation.get_instance;
import diffusion_simulation.many_rand_samples_diff;
import diffusion_simulation.compute_time_instances;

% Set the diffusion jumping distance equal to the pixel size
jumpDist = pixelSize;

% Strategy: simulate intial fragment up until time for first cut, then
% simulate the diffusion for the resulting two fragments up until time for
% second cut etc... 

% extract all lengths of all fragments at each time interval and compute
% the real length in nanometers

fragLengths = extract_frag_lengths(cutInfo.fragmentCollection,intialLength,salt);


% Intiate loop for simulating each time interval


if cutInfo.cutTimes == 0
    cutTimes = maxTime;
else
    cutTimes = [cutInfo.cutTimes,maxTime];
end


% Compute all time instance of interest with the correct resolution
timeInstances = compute_time_instances(timeResolution,maxTime);

% For each instance we need to store the position of the fragments
allFragPos = cell(1,length(timeInstances));
allFragLengths = cell(1,length(timeInstances));
% Set counter for the instance we are situated in
instCount = 1;

% Set time counter
t = 0;

% Set a large number for batch generation of random sampling 
largeNumber = 1000000;
% Count variable to iterate through batch sequence 
step = largeNumber;




for i=1:length(cutTimes)
    
    % calculate the rates for jumping left and right respectivly for each
    % fragment    
    [hopRatesLeft,hopRatesRight] = get_rates(fragLengths{i},diffusionConst,fragLengths{1},pixelSize);    
    mu = sum(hopRatesLeft) + sum(hopRatesRight);
    
    % set maximum simulation time up until next cut
    tMax = cutTimes(i);
    
    % Always generate new random sample when cut happened 
    step = largeNumber;
    
    % Assign positions of all the fragments after cut
    
    if i == 1  % Special case if we are at intial fragment
        
        % create time and postion vector to store evolution in
        fragPos = 0;
        currFragLengths = fragLengths{i};
        
    else % If not in first instance we need to calculate new fragment positions 
        
        
        % Save previous last positions
        prevPos = fragPos;
        
        % Extract previous last lengths and current lenghts
        prevLengths = fragLengths{i-1};
        currFragLengths = fragLengths{i};
        
        % Extract the number for the cutted fragment
        cuttedFragNr = cutInfo.cuttedFragmentNr(i-1);
        
        % Determine new positions
        fragPos = assign_new_pos(prevLengths,currFragLengths,prevPos,cuttedFragNr);
        
    end
    
    
    
    
    while t<tMax
        
        % Generate new random samples give current mu and hop rates if
        % needed 
        
        if step == largeNumber 
            
           [tauVec,fragToMoveVec,leftOrRightVec] = many_rand_samples_diff(...
               hopRatesLeft,mu,largeNumber);
           step = 1;
            
        end
        
        % Draw a time until next move using gillespie
        t = t + tauVec(step);
        
        % Sample the fragment to be moved
        theFragToMove = fragToMoveVec(step);
        
        % Sample if the fragment should be moved left or right
        leftOrRight = leftOrRightVec(step);   
        
        % Used current random number so increment to use next set
        step = step + 1;
        
        
        % Save previous positions
        prevPos =  fragPos;
        
        % Move the sampled fragment in sampled direction if possible
        fragPos = move_frag_if_possible(prevPos,currFragLengths,...
            leftOrRight,theFragToMove,jumpDist); 
        
        % Check if we should save this instance
       
        if t >= timeInstances(instCount)
            
            % Get which time instance we are at (possibly passed multiple)            
            instanceIndex = get_instance(t,timeInstances,instCount);
            
            % If we only passed one instance
            if instanceIndex - instCount < 1
                
                allFragPos{instCount} = fragPos;
                allFragLengths{instCount} = currFragLengths;
                
            else % If we passed multiple instances 
                
                for inst=instCount:1:instanceIndex
                    
                    allFragPos{inst} = fragPos;
                    allFragLengths{inst} = currFragLengths;
                    
                end
                
            end
            
            % Increment the instance counter
            instCount = instanceIndex + 1;

        end        
        
    end
    
end




end

