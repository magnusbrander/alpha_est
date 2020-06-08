function [timeSeriesObs,timeSeriesTrue,nrVisCuts,timeInstances,obsCutTimes, trueCutTimes]...
    = generate_sim_timeseries_lim_time_res(alpha,constVal)

% #########################################################################
%
% FUNCTION: [timeSeriesObs,timeSeriesTrue] = generate_sim_timeseries_lim_time_res(...
%    alpha,constVal)
%
% Description: 
%
% Input: alpha = nicking rate per length, constVal = struct with experiment
% specific constants that does not change when alpha change 
%
% Output: timeSeriesObs = cumulative max for observed time series, 
% timeSeriesTrue = cumulative max for true time series, timeInstances = all
% times in discretized format with spacig of give time resolution up until
% the maximum time, obs/rueCutTimes = the observed/true times when all cuts
% happned
%
% #########################################################################



% Import needed functions
import nicking_simulation.nicking_process_time;
import diffusion_simulation.diffusion_process_limited_res_time;
import diffusion_simulation.plot_diffusion_process;
import diffusion_simulation.merge_fragments_lim_res;
import timeseries_folder.extract_timeseries_lim_time_res;


% Extract constant values 
L = constVal.L; %[basepairs]
pixelSize = constVal.pixelSize; %[nanometers]
frayingDist = constVal.frayingDist; %[basepairs]
diffusionConst = constVal.diffConst; %[nanometers^2 /s]
salt = constVal.salt;
n0 = 2*L*constVal.intialNickDensity;
timeResolution = constVal.timeResolution;
maxTime = constVal.maxTime;


% Compute resolution limit
waveLength = 509; %[nanometers]
NA = 1.4;
sigmaPSF = 1.22*waveLength/(2*NA);
resolvLim = 2*sigmaPSF;


% Perform nicking process of DNA up until given maximum time 
cutInfo = nicking_process_time(L,frayingDist,maxTime,alpha,n0);


% Perform diffusion simulation with limited time resolution up until given time
[allFragPos,allFragLengths,timeInstances]=diffusion_process_limited_res_time(...
    cutInfo,maxTime,timeResolution,L,pixelSize,diffusionConst,salt);

% Merge all fragments below the resolution limit and exclude too short ones
[allFragPosMerged,~ ] = merge_fragments_lim_res(allFragPos,...
    allFragLengths,resolvLim);

% Get the observed and true cumulative max time series 
[timeSeriesObs,timeSeriesTrue,obsCutTimes,nrVisCuts] = ...
    extract_timeseries_lim_time_res(allFragPosMerged,cutInfo,timeInstances);


% Get true times for cuts
trueCutTimes = cutInfo.cutTimes;



end

