function hScoreAndLengthDiff = hScore_and_lengthDiff(nrElement,gradMatrix,scoreMatrix,barcodeWithBinaryEdges,wd,type)



% Function to exract the edge score value and difference in position
% (length) for a pair of edges.

% The function can compute the score and length difference for all posible
% pairs and the nerest neighbour edges, only.



switch type
    
    
    case 'all'
        
        hScoreAndLengthDiff = zeros(nrElement,2);
        
        count = 1;
        
        for i=1:length(gradMatrix)
            
            tempGradMat = gradMatrix{i};
            tempScoreMat = scoreMatrix{i};
            
            
            for time=1:size(tempGradMat,1)
                
                negPositions = find(tempGradMat(time,:) == -1);
                posPositions = find(tempGradMat(time,:) == 1);
                
 
                    if ~isempty(negPositions) && ~isempty(posPositions)

                        negScore = tempScoreMat(time,negPositions);
                        posScore = tempScoreMat(time,posPositions);
                        
                        positionCompare = repmat(posPositions,length(negPositions),1).' - ...
                            repmat(negPositions,length(posPositions),1);
                        
                        scoreCompare = repmat(posScore,length(negScore),1).' + ...
                            repmat(negScore,length(posScore),1);
                        
                        positionCompare = positionCompare(:);
                        scoreCompare = scoreCompare(:);
                        
                        hScoreAndLengthDiff(count:(length(scoreCompare)+count-1),1) = scoreCompare;
                        hScoreAndLengthDiff(count:(length(positionCompare)+count-1),2) = positionCompare;
                        
                        count = count + length(positionCompare);
                        
                    end
                                
            end
            
        end
        
        hScoreAndLengthDiff = hScoreAndLengthDiff(1:count-1,:);
        
    case 'training'
        
        % In the traning case we need to take into acount that it
        % does not have to be a neigbouring edge that we are looking at
        % since some edges were excluded earlier in the data
        
        hScoreAndLengthDiff = zeros(nrElement,2);
        
        count = 1;
        
        for i=1:length(gradMatrix)
            
            tempGradMat = gradMatrix{i};
            tempScoreMat = scoreMatrix{i};
            tempEdgeBinaryImage = barcodeWithBinaryEdges{i};
            
            
            for time=1:size(tempGradMat,1)
                
                
                edgesInBar = find(tempEdgeBinaryImage(time,:) ~= 0);
                
                if length(edgesInBar) == 2
                    leftEdge = edgesInBar(1) - wd;
                    rightEdge = edgesInBar(2) + wd;
                end
                
                % Use the edges to divide the zero-crossings into one left
                % and right set
                
                negPositions = find(tempGradMat(time,:) == -1);
                posPositions = find(tempGradMat(time,:) == 1);
                
                if ~isempty(negPositions) && ~isempty(posPositions)
                    
                        leftNegPositions = negPositions(negPositions<leftEdge);
                        leftPosPositions = posPositions(posPositions<leftEdge);
                        
                        rightNegPositions = negPositions(negPositions>rightEdge);
                        rightPosPositions = posPositions(posPositions>rightEdge);
                        
                        
                        % We need to start with a neg zero cross so therefore we
                        % exclude any positive zero-cross to the left of the first
                        % negative one
                        
                        
                        if ~isempty(leftNegPositions) && ~isempty(leftPosPositions)
                            
                            leftPosPositions = leftPosPositions(leftPosPositions>leftNegPositions(1));
                            
                            if ~isempty(leftPosPositions)
                                for index=1:min(length(leftNegPositions),length(leftPosPositions))
                                    
                                    scoreCompare = tempScoreMat(time,leftNegPositions(index)) +...
                                        tempScoreMat(time,leftPosPositions(index));
                                    positionCompare =  leftPosPositions(index) - leftNegPositions(index);
                                    
                                    hScoreAndLengthDiff(count,1) = scoreCompare;
                                    hScoreAndLengthDiff(count,2) = positionCompare;
                                    
                                    count = count + 1;
                                    
                                    
                                end
                            end
                            
                        end
                        
                        
                        if ~isempty(rightNegPositions) && ~isempty(rightPosPositions)
                            rightPosPositions = rightPosPositions(rightPosPositions>rightNegPositions(1));
                            if ~isempty(rightPosPositions)
                                
                                for index=1:min(length(rightNegPositions),length(rightPosPositions))
                                    
                                    scoreCompare = tempScoreMat(time,rightNegPositions(index)) +...
                                        tempScoreMat(time,rightPosPositions(index));
                                    positionCompare =  rightPosPositions(index) - rightNegPositions(index);
                                    
                                    hScoreAndLengthDiff(count,1) = scoreCompare;
                                    hScoreAndLengthDiff(count,2) = positionCompare;
                                    
                                    count = count + 1;
                                    
                                end
                                
                            end
                            
                        end
                end
                
            end
            
        end
        
        
        
        hScoreAndLengthDiff = hScoreAndLengthDiff(1:count-1,:);
        
        
        
        
    case 'neighbouring'
        
        % In the neighbouring case we have to pair zero-cross edges
        % starting from a negative zero-crossing
        
        
        hScoreAndLengthDiff = zeros(nrElement,2);
        
        count = 1;
        
        for i=1:length(gradMatrix)
            
            
            tempGradMat = gradMatrix{i};
            tempScoreMat = scoreMatrix{i};
            
            for time=1:size(tempGradMat,1)
                
                negPositions = find(tempGradMat(time,:) == -1);
                posPositions = find(tempGradMat(time,:) == 1);
                
                if ~isempty(negPositions) && ~isempty(posPositions)
                    
                    posPositions = posPositions(posPositions>negPositions(1));
                    
                    if ~isempty(posPositions)
                        for index=1:min(length(negPositions),length(posPositions))
                            
                            scoreCompare = tempScoreMat(time,negPositions(index)) +...
                                tempScoreMat(time,posPositions(index));
                            positionCompare =  posPositions(index) - negPositions(index);
                            
                            hScoreAndLengthDiff(count,1) = scoreCompare;
                            hScoreAndLengthDiff(count,2) = positionCompare;
                            
                            count = count + 1;
                            
                            
                        end
                        
                    end
                    
                    
                end
                
            end
            
        end
        
        hScoreAndLengthDiff = hScoreAndLengthDiff(1:count-1,:);
       
        
end





end
