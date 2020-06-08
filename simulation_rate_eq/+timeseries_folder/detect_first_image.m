function [startFrame,framesShort] = detect_first_image(frames)

% Function that detects the first frame of interest in a series of 
% images where the first few frames are dark. 

fprintf('\nDetecting first image of interest...\n')
for i = 1:length(frames)-1
    firstFrame = mean2(frames{i,1}(:))+3*std2(frames{i,1}(:));
    secondFrame = mean2(frames{i+1,1}(:));
     if secondFrame > firstFrame
         startFrame = i+1;
         break
     end
     
end

fprintf('Starting at frame %i.\n',startFrame);

framesShort = cell([],1);
for i = 1:length(frames)-startFrame
    framesShort{i,1} = frames{i+startFrame,1};
end


end
