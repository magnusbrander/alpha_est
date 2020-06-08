function sigIm = add_photons(prevPos,prevTime,currentTime,currFragLengths,...
    timeResolution,boxSizeRow,boxSizeCol,pixelSize,lambda)

% ########################################################################
%
% FUNCTION: output = add_photons(input) 
%
% Description:
%
% Input: 
%
% Output: 
%
% ########################################################################


% Import needed functions
import image_generate.gen_photon_image;

% Compute exposure time up until current move
timeDiff = currentTime - prevTime;
timeFrac = timeDiff/timeResolution;

% Compute previous positions in terms of pixels from left edge of image
prevPosObs = round(prevPos/pixelSize) + round(boxSizeCol/2);

% Get total number of signal pixels
nrSigPix = round(currFragLengths/pixelSize);
totNrSigPix = sum(nrSigPix);


% Create position vectors for all signal pixels 
sigPixCol = zeros(1,totNrSigPix);

% Count variable for all signal pixels
pixCount = 0;

% Loop through all fragments to place signal pixels at right rows 
for frag=1:length(nrSigPix)
    
    % Get ranges for fragment edges
    minCol = prevPosObs(frag) - round(nrSigPix(frag)/2);
    maxCol = prevPosObs(frag) + round(nrSigPix(frag)/2);
    
    % If fragment is inside of box then proceed to place signal pixels
    if minCol > 1 && maxCol < boxSizeCol+1
        
        for p=minCol:1:maxCol
            
            pixCount = pixCount + 1;
            sigPixCol(pixCount) = p;
            
            
        end        
    end   
    
end


% clean up signal pixel array if some fragments were located outside 
sigPixCol = sigPixCol(1:pixCount);
sigPixRow = round(boxSizeRow/2)*ones(1,length(sigPixCol));



sigIm = gen_photon_image(boxSizeRow,boxSizeCol,sigPixRow,sigPixCol,timeFrac,lambda);



end