function output = image_segmentation_pragmatic(image, sigma, userLowLim)

import segementation_util.walk_dist_calc;
import segementation_util.edge_score;
import segementation_util.mol_filt;

% Filter image with LoG filter
n = ceil(6*sigma);
n = n + 1 -mod(n,2);
filt = fspecial('log',n,sigma);
if ~isa(image,'double')
	image = double(image);
end


logim = imfilter(image,filt);


% Find zero crossing contours
thedges = imbinarize(logim,0);
[B,L] = bwboundaries(thedges,'holes');


% Calculate edge score around all edges
[~,Gdir] = imgradient(logim);
dist = walk_dist_calc(sigma); % Very simple function - perhaps do elsewhere.
meh = zeros(1,length(B));
stat = @(h) mean(h);  % This should perhaps be given from the outside


for k = 1:length(B) % Filter out any regions with artifacts in them
    
    meh(k) = edge_score(B{k},logim,Gdir,dist,stat)/sigma^3;
    
end


%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%
lowLim = exp(userLowLim);   % Set the low score threshold to consider a region "signal" (very important)
highLim = exp(24);          % Arbitrary higher bound not utilized at this point. (ignore this for now)
elim = .7;                  % Set lower limit for eccentricity of region (removes dots and circles and keeps long shapes)
ratlim = .3;                % Set lower limit for ratio of area of region to the convex region formed around (removes "wiggly" regions)
lengthLims = [5 1000]; 
%%%%%%%%%%%%%%% TWEAK PARAMETERS %%%%%%%%%%%%%%%



% Filter molecules via the tweak parameters above
accepted = 0;
D = B;
newL = zeros(size(L));
trueedge = cell(1,length(B));

for k = 1:length(B) % Filter any edges with lower scores than lim
	acc = mol_filt(B{k},meh(k),lowLim,highLim,elim,ratlim,lengthLims);
	if acc
		accepted = accepted + 1;
	    trueedge{accepted} = D{k};
		newL(L == k) = accepted;
	else		
		D{k} = [];
	end
end


D = D(~cellfun('isempty',D)); % Remove empty entries (where molecules have been filtered)



output.molecules = D;
output.meh = meh; 


end
