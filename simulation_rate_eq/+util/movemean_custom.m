function outputVec = movemean_custom(inputVec,nrM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% outputVec = movmean_custom(inputVec,nrM)
% 
%
% Description: 
%
% Input: 
%
% Output:  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




outputVec = 0*inputVec;

halfNrM = round(nrM*0.5)-1;

for i=1:length(inputVec)
    
    startPoint = i-halfNrM;
    stopPoint = i+halfNrM;
    
    count = 0;
    
    for j=startPoint:1:stopPoint
        
        if j>0 && j<=length(inputVec)
            outputVec(i) = outputVec(i) + inputVec(j);
            count = count + 1;
        end
        
    end
    
    outputVec(i) = outputVec(i)/count;
 
end

