
function [score,matchIndex] = intact_score(fragToTrack,image,time)



% Extract position for the fragment of interest in last
% instance 
prevPos = fragToTrack(size(fragToTrack,1),1);

% Estimate the length of the fragment 
fragLenEst = sum(fragToTrack(:,2))/length(fragToTrack(:,2));


% Extract current postions and lengths

currPositions = find(image(time,:)~=0);
currLengths = image(time,currPositions);
orginalLength = image(1,image(1,:)~= 0);

% Check if there is a fragment in same position and return a logical value
% as well as which fragments matched the position



import timeseries_folder.position_score
[posScore,posIndMatch] = position_score(prevPos,currPositions,fragLenEst,orginalLength);



% Check if the fragments in the same postion match in length by returning
% a logical array for posIndMatch indicating matches (0 or 1)

if posScore == 1
    
    import timeseries_folder.length_score
    lenIndexScore = length_score(fragLenEst,currLengths,posIndMatch);
    
    % Find the fragments that matched in both length and position
    bothMatch = find(lenIndexScore~=0);
    
    % Extract the specific fragment number that matched
    
    if sum(bothMatch) ~= 0
        
        fragNumThatMatch = posIndMatch(bothMatch);
        
        % Return the index for the matching fragment
        matchIndex = fragNumThatMatch(1);
        
        % Set score to positive 
        score = 1;
        
    else
        score = 0;
        matchIndex = [];
    end
    
    
    
else
    score = 0;
    matchIndex = [];
    
end

    
end

