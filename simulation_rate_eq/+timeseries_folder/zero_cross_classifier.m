function [hScoreMatrix,gradValMatrix] = zero_cross_classifier(logBarcode,zeroCrossLoGBar,barcodeWithBinaryEdges,wd,maxTime,excludeEdges)


inEdgeNotifyCount = 0;
startCol = wd +1;
stopCol = (size(logBarcode,2)-wd-1);
hScoreMatrix = NaN(maxTime,size(logBarcode,2));
gradValMatrix = zeros(maxTime,size(logBarcode,2));



% Chose to include or exclude edges in the computing of edge score
switch excludeEdges
    
    case 'exclude'
        edgeTolerance = wd;        
    case 'include'
        edgeTolerance = 0;
    otherwise
        msg = 'Error occurred.';
        error(msg)
end


for time=1:maxTime    
    
    edgesInBar = find(barcodeWithBinaryEdges(time,:) ~= 0);
    if length(edgesInBar) == 2
        leftEdge = edgesInBar(1);
        rightEdge = edgesInBar(2);       
    end
    
    prevBinVal = zeroCrossLoGBar(time,wd-1);
    
    for col = startCol:stopCol
        
        currentBinVal = zeroCrossLoGBar(time,col);
        
        if abs(col - leftEdge) < edgeTolerance || abs(col - rightEdge) < edgeTolerance   % Estimate if we are located within edge tollerance
            inEdgeNotify = 1;            
        else
            inEdgeNotify = 0;
        end
        
        
        if (currentBinVal-prevBinVal) == 1 && inEdgeNotify == 0
            % This means that we went from negative to positive
            
            % Addjust walking for negative -> positive zero-cross
            tempColNeg =(col-wd):(col-1); % Negative side
            tempColPos = col:(col+wd-1); % Positive side
            
            % Extract the sum in both directions
            negVal = logBarcode(time,tempColNeg);
            posVal = logBarcode(time,tempColPos);
            
            % Check if we crossed another zero-cross while walking
            posValOnNegSide = find(negVal>0,1,'last'); % Take the positive value closest to the zero cross
            negValOnPosSide = find(posVal<0,1,'first'); % Take the negative value closest to the zero cross
            
            % If that is the case limit the sum
            if ~isempty(posValOnNegSide) || ~isempty(negValOnPosSide)
                if ~isempty(posValOnNegSide) && ~isempty(negValOnPosSide)
                    
                    if (wd-posValOnNegSide) > (negValOnPosSide) 
                        % Limit the sum according to the positive side
                        posVal = posVal(1:negValOnPosSide-1);
                        negVal = negVal(wd-length(posVal)+1:end);
                        
                    else
                        % Limit the sum according the negative side
                        negVal = negVal(posValOnNegSide+1:end);
                        posVal = posVal(1:length(negVal));
                       
                        
                    end
                    
                elseif ~isempty(posValOnNegSide) && isempty(negValOnPosSide)
                    % Limit the sum acording the negative side
                     negVal = negVal(posValOnNegSide+1:end);
                     posVal = posVal(1:length(negVal));
                    
                else 
                    % Limit the sum on the postive side
                    posVal = posVal(1:negValOnPosSide-1);
                    negVal = negVal(wd-length(posVal)+1:end);
                    
                end                   
            end
            
    
            
            % Compute the sums to get a score
            negSum = sum(negVal);
            posSum = sum(posVal);
            hScoreMatrix(time,col) = posSum-negSum;
            gradValMatrix(time,col) = currentBinVal-prevBinVal;
            prevBinVal = currentBinVal;
            
        elseif (currentBinVal-prevBinVal) == -1 && inEdgeNotify == 0
            % This means that we went from positive to negative
            
            
            % Addjust walking for positve -> negative zero-cross
            tempColPos =(col-wd):(col-1); 
            tempColNeg = col:(col+wd-1);
            
            % Extract the sum in both directions
            negVal = logBarcode(time,tempColNeg);
            posVal = logBarcode(time,tempColPos);
            
            % Check if we crossed another zero-cross while walking
            posValOnNegSide = find(negVal>0,1,'first');
            negValOnPosSide = find(posVal<0,1,'last');
            
            % If that is the case limit the sum
            if ~isempty(posValOnNegSide) || ~isempty(negValOnPosSide)
                
                if ~isempty(posValOnNegSide) && ~isempty(negValOnPosSide)
                    
                    if (posValOnNegSide) < (wd-negValOnPosSide)
                        % Limit the sum acording to the negative side
                        negVal = negVal(1:posValOnNegSide-1);
                        posVal = posVal(wd-length(negVal)+1:end);
                        
                    else
                        % Limit the sum acording the positive side
                        posVal = posVal(negValOnPosSide+1:end);
                        negVal = negVal(1:length(posVal));
                        
                    end
                    
                elseif ~isempty(posValOnNegSide) && isempty(negValOnPosSide)
                    % Limit the sum acording to the negative side
                     negVal = negVal(1:posValOnNegSide-1);
                     posVal = posVal(wd-length(negVal)+1:end);
                    
                else
                    % Limit the sum on the postive side
                    posVal = posVal(negValOnPosSide+1:end);
                    negVal = negVal(1:length(posVal));
                    
                end
            end
               
            
            % Compute the sums to get a score
            negSum = sum(negVal);
            posSum = sum(posVal);
            hScoreMatrix(time,col) = posSum-negSum;
            gradValMatrix(time,col) = currentBinVal-prevBinVal;
            prevBinVal = currentBinVal;
            
            
            
        elseif abs(prevBinVal-currentBinVal) == 1 && inEdgeNotify == 1
            inEdgeNotifyCount = inEdgeNotifyCount + 1;
            prevBinVal = currentBinVal;           
            
        else
            % This means that we did not cross a zero crossing
            prevBinVal = currentBinVal;
            
        end
            
    end
    
   
end


end
