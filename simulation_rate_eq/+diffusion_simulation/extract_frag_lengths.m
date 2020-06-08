function fragLengths = extract_frag_lengths(fragmentCollection,L,salt)

% FUNCTION: fragLengths = extract_frag_lengths(fragmentCollection,L)
% 
% This funtion extracts the lengths of all fragments in all time intervalls
% and inserts them into a convenient form of cell array ordered such that
% they appear in the spatially correct way in the array
% "fragLengths" has the form of one cell for each time intervall and each
% cell gives the lengths of each fragment in a normal array


% We define the constant Q = 0.9338*(l_K*w_eff/(D_w*D_H))^(1/3), which
% depends upon salt concentration 

switch salt 
    
    case 0.05
        
        Q = 0.9338*((154.4*26)/(100*150))^(1/3);
        
    case 0.5
        
        Q = 0.9338*((116.5*10)/(100*150))^(1/3);
        
    case 2
        
        Q = 0.9338*((105.9*6.2)/(100*150))^(1/3);
        
    case 5
        
        Q = 0.9338*((101.2*4.6)/(100*150))^(1/3);
        
end
% fprintf('---------------')
% isempty(fragmentCollection{1})
% fragmentCollection{:}
% Take into account the extension of the contour length due to dye load
% with approximately 0.44nm per dye molecule (assuming 1:5 staining here)
dyeRatio = 1/5;
dyeLength = 0.44;

% conversion fractor from basepair length to nanometer
C = 0.35;


if isempty(fragmentCollection{1})
    
    fragLengths = cell(1,1);
    fragLengths{1} = Q * L*(C + dyeRatio*dyeLength);
    
else
    
    fragLengths = cell(1,length(fragmentCollection)+1);
    fragLengths{1} = Q * L*(C + dyeRatio*dyeLength);
    
    
    for i=2:length(fragLengths)
        
        tempL = zeros(1,i);
        
        for j=1:i
            
            tempL(j) = size(fragmentCollection{i-1}{j},1);
            tempL(j) = Q * tempL(j)*(C + dyeRatio*dyeLength); % Compute observed length
        end
        
        fragLengths{i} = tempL;
        
    end
    
    
end




end







