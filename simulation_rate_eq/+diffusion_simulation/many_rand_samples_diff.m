function [tauVec,fragToMoveVec,leftOrRight] = many_rand_samples_diff(hopRates,mu,nrSamples)


% ########################################################################
%
% FUNCTION: output = many_rand_samples_diff(hopRates,mu,nrSamples)
%
% Description: 
%
% Input: hopRates = the current hop rates for all fragments in one given 
% direction, mu = the sum of all hop rates in both left and right
% direction, nrSamples = a large number of samples used to speed up the 
% computation of output quantities 
%
% Output: tauVec = waiting time vector, fragToMoveVec = vector with sampled 
% fragments to move , leftOrRight = vector with samples for direction 
%
% ########################################################################



% Compute nrSamples of randomly sampled waiting times given the current mu
% value 
tauVec = - log(rand(1,nrSamples))/mu;

% Compute nrSamples of randomly sampled fragment numbers given the hop
% rates 
fragToMoveVec = randsample(length(hopRates),nrSamples,'true',hopRates);

% Sample nrSamples left or right jumping directions
leftOrRight = randi(2,1,nrSamples);



end
