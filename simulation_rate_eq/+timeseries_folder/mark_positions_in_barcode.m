function barcodeWithPosition = mark_positions_in_barcode(barcode,moleculeInBigBoxesIndex,boxBoundariesForMolecules,bigBoxes,type)

% This function marks postions in a barcode, such as edges or mass
% center
%
% Input:
%   barcode = the barcode itself either in grayscale or rgb format
%
%   moleculeInBigBoxesIndex = the index of molecules enclosed by the
%   current big box of interest
%
%   boxBounariesForMolecules = the boundaries for all detected molecules in
%   the current frame
%
%   bigBoxes = the position of the big box of interest
%
%   type = argument which decides the type of marking: 'binaryEdges', 'centerMass'
%          , 'colourCenterMass', 'colourEdges'
% Ouput:
%   barcodeWithPosition = a barcode marked with the wanted
%   positions



switch type
    
    % Place binary edges on top of black image
    case 'binaryEdges'
        
        barcodeWithPosition = barcode*0;
        
        lengthOfBarcode = size(barcode,1);
        colBigBox = [bigBoxes(1),bigBoxes(1)+bigBoxes(3)];
        
        for j=1:lengthOfBarcode
            
            
            % Extracted molecules of interest for the current frame in the
            % current big box
            enclosedMolecules = moleculeInBigBoxesIndex{j};            
            smallBoxes = boxBoundariesForMolecules{j};
            
            if ~isempty(enclosedMolecules)  
                
                % Extract columns for the smaller enclosed molecules
                columnSmallBoxes = [smallBoxes(enclosedMolecules,1),...
                    smallBoxes(enclosedMolecules,1)+smallBoxes(enclosedMolecules,3)];
                
                col = [columnSmallBoxes(:,1)-colBigBox(1), columnSmallBoxes(:,2) - colBigBox(1)] ;
                
                col = col(:)+1; % Correct for position-pixel discrepancy (zero difference means pixel 1 etc..)
                
                col = sort(col);
                
                % Assign all edge pixels a value of 1 in the barcode 
                barcodeWithPosition(j,col) = 1;                
            end
 
        end
        

        
        
        % Place binary center of mass positions on a black image
    case 'centerMass'
        
        barcodeWithPosition = barcode*0;
        
        lengthOfBarcode = size(barcode,1);
        colBigBox = [bigBoxes(1),bigBoxes(1)+bigBoxes(3)];
        
        for j=1:lengthOfBarcode
            
            
            % Extracted molecules of interest for the current frame in the
            % current big box
            enclosedMolecules = moleculeInBigBoxesIndex{j};            
            smallBoxes = boxBoundariesForMolecules{j};
            
            if ~isempty(enclosedMolecules)  
                
                % Extract columns for the smaller enclosed molecules
                columnSmallBoxes = [smallBoxes(enclosedMolecules,1),...
                    smallBoxes(enclosedMolecules,1)+smallBoxes(enclosedMolecules,3)];
                
                col = [columnSmallBoxes(:,1)-colBigBox(1), columnSmallBoxes(:,2) - colBigBox(1)] ;
                
                col = col(:)+1; % Correct for position-pixel discrepancy (zero difference means pixel 1 etc..)
                
                col = sort(col);
                
                % Compute center of mass column values
                centerCol = round(0.5 * (col(1:2:end)+col(2:2:end)));
                
                % Insert center of mass postions into barcode
                barcodeWithPosition(j,centerCol) = 1; 
                
            end
            
        end
        
        
        
       % Place a coloured pixel at the center of mass on the original
       % barcode
    case 'colourCenterMass'
        
        % Convert barcode to rgb format 
        import timeseries_folder.gray2rgb;
        barcodeWithPosition = gray2rgb(barcode);        
        lengthOfBarcode = size(barcode,1);
        colBigBox = [bigBoxes(1),bigBoxes(1)+bigBoxes(3)];
        
        for j=1:lengthOfBarcode
            
            
            % Extracted molecules of interest for the current frame in the
            % current big box
            enclosedMolecules = moleculeInBigBoxesIndex{j};            
            smallBoxes = boxBoundariesForMolecules{j};
            
            if ~isempty(enclosedMolecules)  
                
                % Extract columns for the smaller enclosed molecules
                columnSmallBoxes = [smallBoxes(enclosedMolecules,1),...
                    smallBoxes(enclosedMolecules,1)+smallBoxes(enclosedMolecules,3)];
                
                col = [columnSmallBoxes(:,1)-colBigBox(1), columnSmallBoxes(:,2) - colBigBox(1)] ;
                
                col = col(:)+1; % Correct for position-pixel discrepancy (zero difference means pixel 1 etc..)
                
                col = sort(col);
                
                % Compute center of mass column values
                centerCol = round(0.5 * (col(1:2:end)+col(2:2:end)));
                
                % Insert center of mass postions into barcode
                for p=1:length(centerCol)
                    
                    barcodeWithPosition(j,centerCol(p),:) = [0,0,1];
                    
                end
                
            end
            
        end
        
 
        %Place coloured edge pixels on the original barcode
    case 'colourEdges'
                
        % Convert barcode to rgb format 
        import timeseries_folder.gray2rgb;
        barcodeWithPosition = gray2rgb(barcode);        
        lengthOfBarcode = size(barcode,1);
        colBigBox = [bigBoxes(1),bigBoxes(1)+bigBoxes(3)];
        
        for j=1:lengthOfBarcode
            
            
            % Extracted molecules of interest for the current frame in the
            % current big box
            enclosedMolecules = moleculeInBigBoxesIndex{j};            
            smallBoxes = boxBoundariesForMolecules{j};
            
            if ~isempty(enclosedMolecules)  
                
                % Extract columns for the smaller enclosed molecules
                columnSmallBoxes = [smallBoxes(enclosedMolecules,1),...
                    smallBoxes(enclosedMolecules,1)+smallBoxes(enclosedMolecules,3)];
                
                col = [columnSmallBoxes(:,1)-colBigBox(1), columnSmallBoxes(:,2) - colBigBox(1)] ;
                
                col = col(:)+1; % Correct for position-pixel discrepancy (zero difference means pixel 1 etc..)
                
                col = sort(col);
                
                % Extract left and right edges
                leftCol = col(1:2:end);
                rightCol = col(2:2:end);
                               
                % Insert edge postions into barcode
                for p=1:length(leftCol)
                    
                    barcodeWithPosition(j,leftCol(p),:) = [0,1,0];
                    barcodeWithPosition(j,rightCol(p),:) = [1,0,0];
                    
                end
                
            end
            
        end
        
        
        
        
        
        
    otherwise
        error('Invalid input mode.')
end


end
