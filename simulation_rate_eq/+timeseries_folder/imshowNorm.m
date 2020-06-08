function imshowNorm(image)



% Create image which can be displed 


minVal = min(min(image));
dispImage = image - minVal;

maxVal = max(max(dispImage));
dispImage = dispImage./maxVal;

figure;
hold on;
imshow(dispImage);


end
