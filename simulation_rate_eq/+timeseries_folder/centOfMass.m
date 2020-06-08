function centerOfMass = centOfMass(positions,weigths)

% This function computes the discrete center of mass for a 1-dimensional
% system 

centerOfMass = (sum(positions.*weigths))/sum(weigths);


end
