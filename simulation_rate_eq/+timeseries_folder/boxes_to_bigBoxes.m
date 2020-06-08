function bigBoxes = boxes_to_bigBoxes(boxes,hExt,vExt)



% A simple function to extend the give boxes by a certain amount 

% hExt is the extention of the boxes in comparison to total length in horisontal (x) direction

% vExt is the extention of the boxes in comparison to total heigth in vertical (y) direction



bigBoxes(:,1) = boxes(:,1) - round((hExt-1)*boxes(:,3)*0.5);
bigBoxes(:,2) = boxes(:,2) - round((vExt-1)*boxes(:,4)*0.5);
bigBoxes(:,3) = round(hExt*boxes(:,3));
bigBoxes(:,4) = round(vExt*boxes(:,4));





end
