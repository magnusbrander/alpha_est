function cutCount = nicking_process_simple(lengthOfDna,frayDist,maxNickNr)


cutCount = zeros(maxNickNr,1);



dna = zeros(lengthOfDna,2);



for nickNr=1:maxNickNr
 
    % Find all the non nicked sites on the DNA 
    [nonNickedSiteRow, nonNickedSiteCol] = find(dna==0);
    
    if isempty(nonNickedSiteRow)
        break;
    end
    
    
    % Choose a random site to nick
    randNumberRow = randi(length(nonNickedSiteRow));
    randSiteRow = nonNickedSiteRow(randNumberRow);
    randSiteCol = nonNickedSiteCol(randNumberRow);
    if randSiteCol == 1
        oppositeCol = 2;
    else
        oppositeCol = 1;
    end
    
    
    % Use 0 for non nickes sites
    % Use 1 for nickes sites
    % Use 2 for excluded zones around cuts
    % Use 3 for cut sites
    
    % Find closest nick on oposite strand 
    nickedSitesOpposite = find(dna(:,oppositeCol) == 1);
    [valDiff,index] = min(abs(nickedSitesOpposite-randSiteRow));
    rowMinDiffPos = nickedSitesOpposite(index);
    % Check if it is within the fraying distance
    if ~isempty(rowMinDiffPos) && valDiff < frayDist
        
        % Place numbers around the cut to indicate taboo zone as below
        
        % ----xxxxxx#xx--------- 
        % ------xx#xxxxxx-------
        % With fraying distance = 4 (x= taboo, # = cut, - = non-nick)
        
        
        % Positive direction opposite strand
        for site=0:frayDist
            
            if (randSiteRow+site) > lengthOfDna || dna(randSiteRow+site,oppositeCol) > 1
                break;
            else 
                dna(randSiteRow+site,oppositeCol) = 2; % Create taboo zone
            end
            
        end

        % Negative direction opposite strand
        for site=1:frayDist
            if (randSiteRow-site) <= 0 || dna(randSiteRow-site,oppositeCol) > 1
                break;
            else
                dna(randSiteRow-site,oppositeCol) = 2; % Create taboo zone
            end
            
        end
        
        
        % Positive direction current strand
        for site=0:frayDist
            
            if (rowMinDiffPos+site) > lengthOfDna || dna(rowMinDiffPos+site,randSiteCol) > 1
                break;
            else
                dna(rowMinDiffPos+site,randSiteCol) = 2; % Create taboo zone
            end
            
        end
        
        
        % Negative direction current strand
        for site=1:frayDist
            
            if (rowMinDiffPos-site) <= 0 || dna(rowMinDiffPos-site,randSiteCol) > 1
                break;
            else
                dna(rowMinDiffPos-site,randSiteCol) = 2; % Create taboo zone
            end  
            
        end

        
        % We now place the cut on the DNA 
        dna(randSiteRow,randSiteCol) = 3;
        dna(rowMinDiffPos,oppositeCol) = 3;
        
        
        % And note that a cut happened at this specific time
        cutCount(nickNr) = 1;
        
    else
        % If no cut can be performed, then place a nick at the selected site 
        dna(randSiteRow,randSiteCol) = 1;
        
    end
    
 
end

















end

