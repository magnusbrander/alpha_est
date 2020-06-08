function moleculeCenterMassAndLengthInBigBox = centOfMassInBigBoxes(moleculeInBigBoxesIndex,boxBoundariesForMolecules,bigBoxes)



moleculeCenterMassAndLengthInBigBox = cell(2,size(moleculeInBigBoxesIndex,2));



for i = 1: size(moleculeInBigBoxesIndex,2)
    
    bigPos = bigBoxes(i,:);
    
   
    for j = 1: size(moleculeInBigBoxesIndex,1) 
        
        
        if ~isempty(moleculeInBigBoxesIndex{j,i})
            
            index = moleculeInBigBoxesIndex{j,i};
            massAndLength = zeros(length(index),2);
            for k=1:length(index)
                
                pos = boxBoundariesForMolecules{j}(index(k),:);
                
                centerMass = round(0.5*pos(3)) + (pos(1)-bigPos(1) + 1);
                
                massAndLength(k,1) = centerMass;
                massAndLength(k,2) = pos(3);
                
            end
            
            moleculeCenterMassAndLengthInBigBox{j,i} = massAndLength;
        end
        
    end
    
end

