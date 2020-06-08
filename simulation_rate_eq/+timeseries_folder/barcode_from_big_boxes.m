function barcodes = barcode_from_big_boxes(frames,bigBoxes,flag)


% This function creates the barcodes for each big box from all frames


barcodes = cell([],1);
[lengthOfBigBoxes,~] = size(bigBoxes);
lengthOfFrames = length(frames);

for i=1:lengthOfBigBoxes
    barcode = zeros(lengthOfFrames,bigBoxes(i,3)); % create barcode image with right dimesions
    for k=1:length(frames)
        tempIm = frames{k};
%         bigBoxes(i,2):(bigBoxes(i,2)+bigBoxes(i,4)-1)
%         bigBoxes(i,1):(bigBoxes(i,1)+bigBoxes(i,3)-1)
        barcode(k,:) = sum(tempIm(bigBoxes(i,2):(bigBoxes(i,2)+bigBoxes(i,4)-1),...
            bigBoxes(i,1):(bigBoxes(i,1)+bigBoxes(i,3)-1)),1); % Compute the mean values in y direction
    end
    % Normalize the barcode 
    
    barcode = barcode./bigBoxes(i,4);
    
    if strcmp(flag,'normalize') 
        barcode = barcode - min(min(barcode));
        barcode = barcode./max(max(barcode));
    end
    
    
    % Save the current barcode
    barcodes{i,1} = barcode;    
    
end






end

