function nrOfElements = nr_of_elements(gradMatrix)



nrOfElements = 0;

for i=1:length(gradMatrix)
    
    tempMat = gradMatrix{i}; 
    
    for j=1:size(tempMat,1)        
        negVal = find(tempMat(j,:) == -1);
        posVal = find(tempMat(j,:) == 1);
        nrOfElements = nrOfElements + length(negVal)*length(posVal);       
    end
    
end


end
