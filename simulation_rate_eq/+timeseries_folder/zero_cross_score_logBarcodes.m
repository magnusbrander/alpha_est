function hScore = zero_cross_score_logBarcodes(logBarcodes,zeroCrossLoGBars,barcodesWithBinaryEdges,wd)


hScore = zeros(1);
inEdgeNotifyCount = 0;
startCol = wd +1;



for i=1:length(zeroCrossLoGBars)
    
    
    
    tempBinLogBar = zeroCrossLoGBars{i};
    tempLogBar = logBarcodes{i};
    stopCol = (size(tempLogBar,2)-wd-1);
    
    edgesInBar = find(barcodesWithBinaryEdges{i}(1,:) ~= 0);
    leftEdge = edgesInBar(1);
    rightEdge = edgesInBar(2);
    
    prevBinVal = tempBinLogBar(1,wd-1);
    
    for col = startCol:stopCol
        
        currentBinVal = tempBinLogBar(1,col);
        
        if abs(col - leftEdge) < 0 || abs(col - rightEdge) < 0   % Estimate if we are located within edge tollerance
            inEdgeNotify = 1;            
        else
            inEdgeNotify = 0;
        end
        
        
        if (currentBinVal-prevBinVal) == 1 && inEdgeNotify == 0
            % This means that we went from negative to positive
            prevBinVal = currentBinVal;
            
            % Addjust walking for negative -> positive zero-cross
            tempColNeg =(col-wd):(col-1); 
            tempColPos = col:(col+wd-1);
            
            % Extract the sum in both directions
            negVal = tempLogBar(1,tempColNeg);
            posVal = tempLogBar(1,tempColPos);
            
            % Check if we crossed another zero-cross while walking
            posValInNeg = find(negVal>0,1,'last');
            negValInPos = find(posVal<0,1,'first');
            
            % If that is the case limit the sum
            if ~isempty(posValInNeg) || ~isempty(negValInPos)
                if ~isempty(posValInNeg) && ~isempty(negValInPos)
                    
                    if (wd-posValInNeg) > (2*wd-negValInPos)
                        % Limit the sum on the negative side
                        negVal = negVal(posValInNeg+1:end);
                        posVal = posVal(1:lenght(negVal));
                        
                    else
                        % Limit the sum on the positive side
                        posVal = posVal(1:negValInPos-1);
                        negVal = negVal(wd-length(posVal):end);
                        
                    end
                    
                elseif ~isempty(posValInNeg) && isempty(negValInPos)
                    % Limit the sum on the negative side
                     negVal = negVal(posValInNeg+1:end);
                     posVal = posVal(1:lenght(negVal));
                    
                else 
                    % Limit the sum on the postive side
                    posVal = posVal(1:negValInPos-1);
                    negVal = negVal(length(negVal)-length(posVal):end);
                    
                end                   
            end
            
    
            
            % Compute the sums to get a score
            negSum = sum(negVal);
            posSum = sum(posVal);
            hScore(end+1) = (posSum-negSum);
            
            
        elseif (currentBinVal-prevBinVal) == -1 && inEdgeNotify == 0
            % This means that we went from positive to negative
            pprevBinVal = currentBinVal;
            
            % Addjust walking for positve -> negative zero-cross
            tempColPos =(col-wd):(col-1); 
            tempColNeg = col:(col+wd-1);
            
            % Extract the sum in both directions
            negVal = tempLogBar(1,tempColNeg);
            posVal = tempLogBar(1,tempColPos);
            
            % Check if we crossed another zero-cross while walking
            posValInNeg = find(negVal>0,1,'first');
            negValInPos = find(posVal<0,1,'last');
            
            % If that is the case limit the sum
            if ~isempty(posValInNeg) || ~isempty(negValInPos)
                
                if ~isempty(posValInNeg) && ~isempty(negValInPos)
                    
                    if (2*wd-posValInNeg) < (wd-negValInPos)
                        % Limit the sum on the negative side
                        negVal = negVal(1:posValInNeg-1);
                        posVal = posVal(wd-length(negVal):end);
                        
                    else
                        % Limit the sum on the positive side
                        posVal = posVal(negValInPos:end);
                        negVal = negVal(1:length(posVal));
                        
                    end
                    
                elseif ~isempty(posValInNeg) && isempty(negValInPos)
                    % Limit the sum on the negative side
                     negVal = negVal(1:posValInNeg-1);
                     posVal = posVal(wd-length(negVal):end);
                    
                else
                    % Limit the sum on the postive side
                    posVal = posVal(negValInPos:end);
                    negVal = negVal(1:length(posVal));
                    
                end
            end
               
            
            % Compute the sums to get a score
            negSum = sum(negVal);
            posSum = sum(posVal);
            hScore(end+1) = (posSum-negSum);
            
            
        elseif abs(prevBinVal-currentBinVal) == 1 && inEdgeNotify == 1
            inEdgeNotifyCount = inEdgeNotifyCount + 1;
            
        else
            % This means that we did not cross a zero crossing
            prevBinVal = currentBinVal;
            
        end
            
    end
    
    
    
end

inEdgeNotifyCount




end
