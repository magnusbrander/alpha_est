function bgIm = gen_bg_image(dim1, dim2,lambda,g)

% ########################################################################
%
% FUNCTION: bgIm = gen_bg_image(dim1, dim2,lambda,g)
%
% Description: Generates a background image with Poisson distributed 
% number of photons in each pixel with mean "lambda". The photons are converted 
% through a EMCCD camera noise model with gain "g" described in 
% "https://doi.org/10.1038/s41592-019-0364-4"
%
% Input: dim1 = number of rows, dim2 = number of columns, lambda = average
% number of photons for all pixels, g = gain 
%
% Output: bgIm = the background image 
%
% ########################################################################

r = 74.4;                  % readout noise
f = 45;				       % photon to electron conversion factor
offset = 100;              % Addition to readout pixel value
c = 0.002;                 % clock-induced charge
QE = 1;                    % quantum efficiency 
maxCount = 65535;          % Maximum pixel count 
sigma_R = (r/f);


% Incoming photons to the pixels 
n_ph = poissrnd(lambda,dim1,dim2);

% Incoming electrons to detector
nie = poissrnd(QE*n_ph + c); 

% Number of output electrons + readout noise
noe = gamrnd(nie,g)+ sigma_R*randn(dim1,dim2);

% Rounding to pixel counts
bgIm = min(floor(noe/f) + offset,maxCount); 


end
