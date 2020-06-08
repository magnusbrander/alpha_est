function gammaIm = gen_gamma_image(photonIm,g)

% Generates random image counts from a Poisson distribution
% transformed through an EMCCD camera. See Methods in 
% Sage, "Super-resolution fight club: assessment of 2D and 3D 
% single-molecule localization microscopy software", Nature Methods (2019).


%g = 254;                           % gain in EMCCD camera
r = 74.4;                          % readout noise
f = 45;                            % analogue-to-digital conversion factor
sigma_R = (r/f); 

maxCount = 65535;
% Incoming electrons to detector

[sigRow,sigCol] = find(photonIm ~= 0); 

noe = gamrnd(photonIm,g);


% Adding readout noise and convert to digital numbers to get
% the image counts (nic) [not concerting to integers, like Sage et al.
% does]

readOutNoise = sigma_R*randn(1,length(sigRow));

for i = 1:length(sigRow)
    
    noe(sigRow(i),sigCol(i)) = noe(sigRow(i),sigCol(i)) + readOutNoise(i);
    
end


gammaIm = floor(noe/f); 

gammaIm = min(gammaIm,maxCount);



end
