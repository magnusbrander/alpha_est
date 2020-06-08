function bgIm = gen_bg_image(dim1, dim2,lambda,g)


r = 74.4;                          % readout noise
f = 45;								% photon to electron conversion factor
offset = 100;						%Addition to readout pixel value
c = 0.002;                         % clock-induced charge
QE = 1;                          % quantum efficiency 
maxCount = 65535;
sigma_R = (r/f); 
                          
% Incoming electrons to detector
nie = poissrnd(QE*lambda + c,dim1,dim2); 


% Number of output electrons 

% Adding readout noise and convert to digital numbers to get
% the image counts (nic) [not concerting to integers, like Sage et al.
% does]

noe = gamrnd(nie,g)+ sigma_R*randn(dim1,dim2);

bgIm = min(floor(noe/f) + offset,maxCount); 


end
