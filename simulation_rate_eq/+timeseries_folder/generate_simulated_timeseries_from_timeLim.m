function [timeSeriesObs,timeSeriesTrue] = generate_simulated_timeseries_from_timeLim(alpha,Xi,L,...
    maxTime,pixelSize,diffusionConst,salt,initialNickNr)


% #########################################################################
%
% FUNCTION: cutTimeDist = generate_cutTime_dist(input)
%
% Description: A funtion for generating, based on simulation, distributions
% for differnce between cut times given some values of nicking rate, 
% fraying distance, length of fragment. 
%
% Input: alpha = nicking rate per length, Xi = fraying distance in
% basepairs, L = length of fragment in basepairs, cutNr = number of cuts,
% pixel size [meters]
%
% Output: timeSeriesObs = the observed time series after acconting for
% noise during detection and limited reslution, timeSeriesTrue = the
% underlying true time series we simulated 
%
% #########################################################################



% --------------- Needed constants ---------------------

% Compute the diffusion jumping length in terms of basepairs 


% Compute resolution limit
waveLength = 509; %[nanometers]
NA = 1.4;
sigmaPSF = 1.22*waveLength/(2*NA);
resolvLim = 3*sigmaPSF ;% + detectionUncertainty_BP;



% -------------- Nicking simulation ------------------ %

% First we need to simulate the nicking process complemented with the
% position of the cuts to allow for recognition of the fragments



import nicking_simulation.nicking_process_time;
[cutInfo,~] = nicking_process_time(L,Xi,maxTime,alpha,initialNickNr);


% -------------- Diffusion simulation ------------------- %

% Run diffusion simulation over the nicking process
import diffusion_simulation.diffusion_process_time;
fragmentEvolution = diffusion_process_time(cutInfo,maxTime,L,pixelSize,diffusionConst,salt);

if isempty(fragmentEvolution)
    timeSeriesObs = 0;
    timeSeriesTrue = 0;
    return;
end

import diffusion_simulation.extract_frag_lengths;
fragLengths = extract_frag_lengths(cutInfo.fragmentCollection,L,salt);



% -------------- Limited resolution --------------------- %

% The idea is to go through all fragments and merge together pairs if they are
% below the resolution limit 

import diffusion_simulation.merge_fragments;
combinedFragments = merge_fragments(fragmentEvolution,resolvLim,fragLengths);


% -------------- Extract time series ------------------- %

% Here we extract both the true and observed time series 

import diffusion_simulation.extract_timeseries_from_simulation;
[timeSeriesObs,timeSeriesTrue] = extract_timeseries_from_simulation(combinedFragments);




end














