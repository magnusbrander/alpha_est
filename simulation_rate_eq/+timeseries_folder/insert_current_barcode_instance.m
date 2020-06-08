function barcodeInstance = insert_current_barcode_instance(image,bigBoxes,tempBoxes,action)


% ###################################################################
%
% FUNCTION: barcodeInstance = insert_current_barcode_instance(image,bigBoxes,tempBoxes,action)
%
% Description: This function generates a single time instance of a barcode
% with different features according to specification
%
% Input: image = current image, bigBoxes = big rectangles constraining region
% of interest for intially detected molecules, tempBoxes = boxes given by
% molecules detected in this time instance, action = specify which type of
% barcode we want to generate
%
% Output: barcodeInstances = the current barcode instances for all intilly
% detected molecules
%
% ###################################################################


barcodeInstance = cell(1,size(bigBoxes,1));


switch action
    
    
    % Create an intensity profile from current image and big boxes
    case 'original'
        
        for i=1:size(bigBoxes,1)
            
            % get position, widht and height of region
            x = bigBoxes(i,1);
            y = bigBoxes(i,2);
            dx = bigBoxes(i,3);
            dy = bigBoxes(i,4);
            tempSmallImage = image(y:(y+dy),x:(x+dx));
            barcodeInstance{i} = sum(tempSmallImage,1);
            
        end
        
        
        % Create a binary profile indicating center of mass for detected
        % molecules
    case 'binaryCM'
        
        % Compute area of all small boxes
        areaOfMol = tempBoxes(:,3).*tempBoxes(:,4);
        % Compute center of mass for all small boxes along x-dimension
        smallBoxesCM = round( 0.5*(2*tempBoxes(:,1)+tempBoxes(:,3)));
        
        for i=1:size(bigBoxes,1)
            
            % Get current big box
            tempBigBox = bigBoxes(i,:);
            binaryProfile = zeros(1,tempBigBox(3)+1);
            
            % Get the normalized small box area enclosed by current big box
            intersection = rectint(tempBoxes,tempBigBox)./areaOfMol;
            
            % Get indices for molecules located within the current big box
            indexMatch = find(intersection==1);
            
            % Insert their center of mass into the binary profile
            if ~isempty(indexMatch)
                for j=1:length(indexMatch)
                    % Position within binary profile for current small
                    % molecule
                    pTemp = smallBoxesCM(indexMatch(j)) - tempBigBox(1) + 1;
                    binaryProfile(pTemp) = 1;
                end
            end
            % Save the current binary profile
            barcodeInstance{i} = binaryProfile;
            
        end
        
        
        
        % Create intensity profile from current image and big boxes with added
        % colour edges marked out for detected molecules
    case 'colourEdges'
        
        % Compute area of all small boxes
        areaOfMol = tempBoxes(:,3).*tempBoxes(:,4);
        
        for i=1:size(bigBoxes,1)
            % Get current big box
            tempBigBox = bigBoxes(i,:);
            
            % get position, width and height of region
            x = tempBigBox(1);
            y = tempBigBox(2);
            dx = tempBigBox(3);
            dy = tempBigBox(4);
            
            % get image corresponding to current big box
            tempSmallImage = image(y:(y+dy),x:(x+dx));
            
            % Compute intensity profile
            intensityProfile = sum(tempSmallImage,1);
            
            % Get the normalized small box area enclosed by current big box
            intersection = rectint(tempBoxes,tempBigBox)./areaOfMol;
            
            % Get indices for molecules located within the current big box
            indexMatch = find(intersection==1);
            
            % Insert their center of mass into the binary profile
            if ~isempty(indexMatch)
                for j=1:length(indexMatch)
                    % Position within binary profile for current small
                    % molecule
                    pTempLeft = tempBoxes(indexMatch(j),1) - tempBigBox(1) + 1;
                    pTempRight = tempBoxes(indexMatch(j),3) + pTempLeft;
                    intensityProfile(pTempLeft) = -inf; % Mark left edges as -inf
                    intensityProfile(pTempRight) = inf; % Mark right edges as inf
                end
            end
            % Save the current binary profile
            barcodeInstance{i} = intensityProfile;
            
        end
        
        
        % Create a binary profile with intensity value given by the length of
        % the fragment placed at the center of mass
    case 'binaryWithLength'
        
        % Compute area of all small boxes
        areaOfMol = tempBoxes(:,3).*tempBoxes(:,4);
        % Compute center of mass for all small boxes along x-dimension
        smallBoxesCM = round( 0.5*(2*tempBoxes(:,1)+tempBoxes(:,3)));
        
        for i=1:size(bigBoxes,1)
            
            % Get current big box
            tempBigBox = bigBoxes(i,:);
            binaryLengthProfile = zeros(1,tempBigBox(3)+1);
            
            % Get the normalized small box area enclosed by current big box
            intersection = rectint(tempBoxes,tempBigBox)./areaOfMol;
            
            % Get indices for molecules located within the current big box
            indexMatch = find(intersection==1);
            
            % Insert their center of mass into the binary profile
            if ~isempty(indexMatch)
                for j=1:length(indexMatch)
                    % Position within binary profile for current small
                    % molecule
                    pTemp = smallBoxesCM(indexMatch(j)) - tempBigBox(1) + 1;
                    LTemp = tempBoxes(indexMatch(j), 3);
                    binaryLengthProfile(pTemp) = LTemp;
                end
            end
            % Save the current binary profile
            barcodeInstance{i} = binaryLengthProfile;
            
        end
        
        
        
        
        
end





end
