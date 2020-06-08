
function [fragments] = fragment_tracking(posLenImage)



fragments = cell(1,1);

im = posLenImage;

% Find the first fragment

pos = find(im(1,:) ~= 0);

len = im(1,pos);

prevState = 1;
fragmentID = 1;
mergedPartnerID = 0;
matchIndex = 1;
time = 1; 

fragments{1,1} = [pos,len,prevState,fragmentID,mergedPartnerID,matchIndex,time];
% Initiate neccerasy array to store scores
nrOfScores = 7;
stateScore = zeros(size(posLenImage,1),nrOfScores);

posInd = 1;
lenInd = 2;
prevStateInd = 3;
fragIDInd = 4;
mergedPartnerInd = 5;
mergeMatchInd = 6;
timeInd = 7;

% Set a threshold value for the total length score under which we stop the
% tracking
epsilon = 0;


for time = 2:size(posLenImage,1)
    
    
    % Here we want to check all possible cases that can happend in order to
    % figure what has happened to the fragment
    
    
    %   1. Check if fragment is intact
    %
    %   2. Check if fragment split but we detected the children
    %
    %   3. Check if fragment is merged with another fragment
    %
    %   4, Check if total length is conserved
    
    % Temporary binary score for the different cases
    tempScore = zeros(1,nrOfScores);
    % Temporary cell array for fragments, used in order to be flexible for
    % merging and split
    tempFragments = cell(1,1);
    % Help variable for fragment number
    fragmentCount = 1;
    % Help variable to know when to add a new row to the fragment cell
    % array
    fragmentSplit = 0;
    
    
    % Current information for the fragments
    currFragPositions = find(im(time,:) ~= 0);
    currFragLengths = im(time,currFragPositions);
    currFragments = [currFragPositions.',currFragLengths.'];
    
    % ############ Import functions ###################
    import timeseries_folder.intact_score
    import timeseries_folder.splitUp_score
    import timeseries_folder.appendNum
    import timeseries_folder.merged_score
    
    tempNumOfFrag = length(fragments(size(fragments,1),:));
    frag = 1;
    while frag <= tempNumOfFrag
        
        
        totFragToTrack = fragments{size(fragments,1),frag};
        fragmentID = totFragToTrack(size(totFragToTrack,1),4);
        prevState = totFragToTrack(size(totFragToTrack,1),3);
        
        
        if prevState == 1
            
            % Check for intact case
            [tempScore(1),matchIndex] = intact_score(totFragToTrack,im,time);
            
            % Check for split case if previous is negative
            if tempScore(1) == 1
                % If the fragment is found intact we simply update the
                % fragment sequence           
                fragUpdate = [currFragPositions(matchIndex),currFragLengths(matchIndex),prevState,fragmentID,0,matchIndex,time];
                totFragToTrack = [totFragToTrack;fragUpdate];
                tempFragments{fragmentCount} = totFragToTrack;
                fragmentCount = fragmentCount + 1;
            else
                % If fragment is not intact check for split up
                [tempScore(2),matchIndex] = splitUp_score(totFragToTrack,im,time);
                
                % If fragment is split up, then proceed by creating new fragments
                if tempScore(2) == 1
                    
                    fprintf('Detcted split-up')
                    time
                    
                    fragmentSplit = 1;
                    
                    % Create the new fragments from the matching pair and
                    % use append number function to indicate the sequence
                    % of parents
                    firstFrag = [currFragPositions(matchIndex(1)),currFragLengths(matchIndex(1)),...
                        prevState,appendNum(fragmentID,1),0,matchIndex(1),time];
                    secondFrag = [currFragPositions(matchIndex(2)),currFragLengths(matchIndex(2)),...
                        prevState,appendNum(fragmentID,2),0,matchIndex(2),time];
                    
                    % Insert the new fragments in the temporary fragment
                    % cell
                    tempFragments{fragmentCount} = firstFrag;
                    fragmentCount = fragmentCount + 1;
                    tempFragments{fragmentCount} = secondFrag;
                    fragmentCount = fragmentCount + 1;
                    
                    
                    
                    
                else
                    % If fragment is not intact neigther split up: check
                    % merged test if we are not checking the last fragment
                    
                    
                    if tempNumOfFrag > frag
                        potenMatchFrag = fragments{size(fragments,1),frag+1};
                        [tempScore(3),matchIndex] = merged_score(totFragToTrack,potenMatchFrag,im,time);
                        
                        % If fragment appears to have merged then
                        if tempScore(3) == 1
                            fprintf('Detected merge')
                            time
                            
                            % Set the previous state to merged state
                            prevState = 3;
                            
                            % We now use the previous and current lengths and the
                            % current position to estimate the true current
                            % position of the merged fragments
                            
                            % Position and length of merged fragment
                            xm = currFragPositions(matchIndex);
                            lm = currFragLengths(matchIndex);
                            
                            % Previous length of the first and second fragment
                            l1 = totFragToTrack(size(totFragToTrack,1),2);
                            l2 = potenMatchFrag(size(potenMatchFrag,1),2);
                            
                            % Estimate the position of the first and second fragment
                            % based on the distance from the left and right
                            % edge, respectively
                            x1 = round(xm - lm*0.5 + l1*0.5);
                            x2 = round(xm + lm*0.5 -l2*0.5);
                            
                            
                            % Updated fragments information
                            fragmentIDFirst = totFragToTrack(size(totFragToTrack,1),4);
                            fragmentIDSecond = potenMatchFrag(size(potenMatchFrag,1),4);
                            
                            fragUpdateFirst = [x1,l1,prevState,fragmentIDFirst,fragmentIDSecond,matchIndex,time];
                            totFragToTrack = [totFragToTrack;fragUpdateFirst];
                            tempFragments{fragmentCount} = totFragToTrack;
                            fragmentCount = fragmentCount + 1;
                            
                            fragUpdateSecond = [x2,l2,prevState,fragmentIDSecond,fragmentIDFirst,matchIndex,time];
                            potenMatchFrag = [potenMatchFrag;fragUpdateSecond];
                            tempFragments{fragmentCount} = potenMatchFrag;
                            fragmentCount = fragmentCount + 1;
                            
                            % No need to check the next fragment now as it
                            % is merged so we increment the loop variable
                            % for the fragment by one extra here
                            frag = frag + 1;
                        else
                            fprintf('Change did not match any of the alternatives')
                            time
                            return;
                            
                        end
                    else
                        fprintf('Situated in last fragment, did not attempt to merge')
                        time
                        return;
                    end
                    
                    
                    
                    
                    
                    
                end
                
                
            end
        elseif prevState == 3
           
            % If we know that the current fragment is located in a merged
            % state we first check if it is still merged. We do this by
            % checking if the merged fragment is intact
            
            % Extract the merged fragment we want to track
            mIndex = totFragToTrack(size(totFragToTrack,1),mergeMatchInd);
            timeOfInterest = totFragToTrack(size(totFragToTrack,1),timeInd);
            prevFragIndex = find(im(timeOfInterest,:)~= 0);
            
            xm = (prevFragIndex(mIndex));
            ml = im(timeOfInterest,xm);
            
            % Create previously detected merged fragment             
            mergedFragToTrack = [xm,ml,0,0,0,0,timeOfInterest];
           
            % Check intact case for merged fragment 
            [tempScore(3),matchIndex] = intact_score(mergedFragToTrack,im,time);
            
            % Check if the merged fragment is still merged 
            if tempScore(3) == 1
                
                fprintf('Constant merged state detected\n')
                % Extract the second fragment in the merged fragment (with 
                % the current assumption it is easy as max two fragments
                % can form a merged fragment so the partner fragment have
                % to be the next fragment in order)                
                partnerMatchFrag = fragments{size(fragments,1),frag+1};
                
                % We now use the previous and current lengths and the
                % current position to estimate the true current
                % position of the merged fragments
                
                
                % Previous length of the first and second fragment
                l1 = totFragToTrack(size(totFragToTrack,1),lenInd);
                l2 = partnerMatchFrag(size(partnerMatchFrag,1),lenInd);
                
                % Estimate the position of the first and second fragment
                % based on the distance from the left and right
                % edge of the merged fragment, respectively
                x1 = round(xm - lm*0.5 + l1*0.5);
                x2 = round(xm + lm*0.5 -l2*0.5);
                
                

                % Updated fragments information
                fragmentIDFirst = totFragToTrack(size(totFragToTrack,1),fragIDInd);
                fragmentIDSecond = partnerMatchFrag(size(partnerMatchFrag,1),fragIDInd);
                
                fragUpdateFirst = [x1,l1,prevState,fragmentIDFirst,fragmentIDSecond,matchIndex,time];
                totFragToTrack = [totFragToTrack;fragUpdateFirst];
                tempFragments{fragmentCount} = totFragToTrack;
                fragmentCount = fragmentCount + 1;
                
                
                fragUpdateSecond = [x2,l2,prevState,fragmentIDSecond,fragmentIDFirst,matchIndex,time];
                partnerMatchFrag = [partnerMatchFrag;fragUpdateSecond];
                tempFragments{fragmentCount} = partnerMatchFrag;
                fragmentCount = fragmentCount + 1;
                
                % No need to check the next fragment now as it
                % is merged so we increment the loop variable
                % for the fragment by one extra here
                frag = frag + 1;
                
            else
                
                % At this point we know that the fragment was previously
                % merged but appears not to be any longer. Therefor, we
                % check it has unmerged and we can again detect both of the
                % individual fragments 
                partnerMatchFrag = fragments{size(fragments,1),frag+1};
                
                % Check for intact case for both fragments
                [scoreFirst,matchIndexFirst] = intact_score(totFragToTrack,im,time);
                % Check for intact case
                [scoreSecond,matchIndexSecond] = intact_score(partnerMatchFrag,im,time);
                
                if scoreFirst == 1 && scoreSecond == 1
                    
                    % Here we know that the fragments were succesfully
                    % detected again so we simply update their positions
                    % and previous state
                    
                    % Current postions and lengths
                    x1 = currFragPositions(matchIndexFirst);
                    x2 = currFragPositions(matchIndexSecond);
                    l1 = currFragLengths(matchIndexFirst);
                    l2 = currFragLengths(matchIndexSecond);
                    
                    % Extract previous IDs for fragments
                    fragmentIDFirst = totFragToTrack(size(totFragToTrack,1),fragIDInd);
                    fragmentIDSecond = partnerMatchFrag(size(partnerMatchFrag,1),fragIDInd);
                    
                    % Set the previous state to non merged 
                    prevState = 1;
                    
                    fragUpdateFirst = [x1,l1,prevState,fragmentIDFirst,0,matchIndexFirst,time];
                    totFragToTrack = [totFragToTrack;fragUpdateFirst];
                    tempFragments{fragmentCount} = totFragToTrack;
                    fragmentCount = fragmentCount + 1;
                    
                    
                    fragUpdateSecond = [x2,l2,prevState,fragmentIDSecond,0,matchIndexSecond,time];
                    partnerMatchFrag = [partnerMatchFrag;fragUpdateSecond];
                    tempFragments{fragmentCount} = partnerMatchFrag;
                    fragmentCount = fragmentCount + 1;
                    
                    % No need to check the next fragment now so we 
                    % increment the loop variable for the fragment by 
                    % one extra here
                    frag = frag + 1;
                    
                    fprintf('Verified unmerging')
                    time
                    
                    
                else
                    fprintf('Did not managed to detected the unmerged fragments')
                    return;
                end
                
            end
            
        else
            fprintf('Invalid state')
            return;
        end
        
        % Increment the loop variable to track the next fragment
        frag = frag + 1;
        
    end
    
    
    
    % Transfer the temporary fragments to the fragments cell array
    if fragmentSplit == 1
        
        row = size(fragments,1) + 1;
        for k=1:length(tempFragments)
            fragments{row,k} = tempFragments{k};
        end
        
    else
        
        row = size(fragments,1);
        for k=1:length(tempFragments)
            fragments{row,k} = tempFragments{k};
        end
        
    end
    
    
    
    stateScore(time,:) = tempScore;
    
    
    if tempScore(4) < epsilon
        return;
    end
    
    
    
end


% Return the fragment sequence


end
