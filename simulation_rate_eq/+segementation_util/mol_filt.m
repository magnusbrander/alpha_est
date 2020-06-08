function acc = mol_filt(B,score,lowLim,highLim,elim,ratlim,lengthLims)
l = sqrt( (max(B(:,2))-min(B(:,2)))^2 + (max(B(:,1))-min(B(:,1)))^2);
lOk = (l > lengthLims(1) && l < lengthLims(2));

import segementation_util.cont_draw
if score > lowLim && score < highLim && lOk
	[~,ecc,aRat] = cont_draw(B);
	testofboundary = (ecc > elim && aRat > ratlim);
	if testofboundary
		acc = true;
	else 
		acc = false;
	end
else
	acc = false;
end
