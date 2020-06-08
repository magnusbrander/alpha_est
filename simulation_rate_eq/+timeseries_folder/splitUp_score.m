function [score,indexMatch] = splitUp_score(fragToTrack,image,time)



% Compare pairs of fragment to the orginal fragment to check if they match
% in position and length

% Determine the average length of the fragment 
avgFragLen = sum(fragToTrack(:,2))/length(fragToTrack(:,2));

prevPos = fragToTrack(size(fragToTrack,1),1);

% Get position of fragments 
currPos = find(image(time,:)~= 0);

currLen= image(time,currPos);


% Get an estimate of intial length
orignalLength = image(1,image(1,:)~= 0);

% Get total number of detected fragments in time frame
nrOfFrag = length(currPos);


score = 0;
indexMatch = [];
% If the number of fragments is equal or larger than two, proceed,
% ohterwise stop the atempt to match fragments
if nrOfFrag >= 2
    
    % ##### Import needed functions ######
    import timeseries_folder.position_score
    import timeseries_folder.length_score
    import timeseries_folder.centOfMass
    % #####################################
    
    i = 1;
    foundMatch = 0;
    while i <= (nrOfFrag-1) && foundMatch == 0
        
        % Compute center of mass position and total length for the pair of
        % fragments 
        commonLen = currLen(i)+currLen(i+1);
        positions = [currPos(i),currPos(i+1)];
        lengths = [currLen(i),currLen(i+1)];
        centerOfMass = centOfMass(positions,lengths);
        % Take the average length of the two fragments as the length
        % estimate in the computation of the diffusion speed
        avgLen = round(mean([currLen(i),currLen(i+1)]));
        
        % Run the position score test on the composite fragment
        [posScore,posIndMatch] = position_score(prevPos,centerOfMass,avgLen,orignalLength);
        
        
        % Check if the length of the composite fragment matches the
        % original fragment if position score is positive  
        
        if posScore == 1
            
            lenIndexScore = length_score(avgFragLen,commonLen,posIndMatch);
            
            if lenIndexScore == 1
                score = 1; 
                indexMatch = [i,i+1];
                foundMatch = 1;
            end
            
            
            
        end
        
       % Increment to move on to the next pair  
       i = i+1;  
        
    end
        
else 
    score = 0;
    indexMatch = [];
    
end

end 