function photonIm = gen_photon_image(nrRowIm,nrColIm, sigRows, sigCols, illuminationFrac,lambda)

% Generates random image counts from a Poisson distribution
% transformed through an EMCCD camera. See Methods in 
% Sage, "Super-resolution fight club: assessment of 2D and 3D 
% single-molecule localization microscopy software", Nature Methods (2019).


% pixelSize = 160;
% waveLength = 509;
% NA = 1.49;
% sigma = 1.22*waveLength/(2*NA) / pixelSize;


noOfRndNumbers = length(sigRows);  % number of random numbers to be generated
%lambda = 300*pi*sigma^2;     % poisson parameter 
                                   % (number of photons, n_photIn in Sage et al.)
c = 0.002;                         % clock-induced charge
QE = 0.9;                          % quantum efficiency 
                          
% Incoming electrons to detector
nie = poissrnd((QE*lambda + c)*illuminationFrac,1,noOfRndNumbers); 

photonIm = zeros(nrRowIm,nrColIm);


for i = 1:length(sigCols)
    
    photonIm(sigRows(i),sigCols(i)) = nie(i);
    
end



end
