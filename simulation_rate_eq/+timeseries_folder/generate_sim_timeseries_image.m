function [timeSeriesObs,timeSeriesTrue,nrVisCuts]= generate_sim_timeseries_image(...
    alpha,constVal)

% #########################################################################
%
% FUNCTION: [timeSeriesObs,timeSeriesTrue] = generate_sim_timeseries_lim_time_res(...
%    alpha,constVal)
%
% Description: 
%
% Input: alpha = nicking rate per length, constVal = struct with experiment
% specific constants that does not change when alpha change 
%0
% Output: timeSeriesObs = cumulative max for observed time series, 
% timeSeriesTrue = cumulative max for true time series, timeInstances = all
% times in discretized format with spacig of give time resolution up until
% the maximum time, obs/rueCutTimes = the observed/true times when all cuts
% happned
%
% #########################################################################



% Import needed functions
import nicking_simulation.nicking_process_time;
import diffusion_simulation.diffusion_process_image_speed;
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
waveLength = constVal.waveLength;
NA = constVal.NA;


% Perform nicking process of DNA up until given maximum time 
cutInfo = nicking_process_time(L,frayingDist,maxTime,alpha,n0);


% Perform diffusion simulation with image segmentation up until given time
[allFragPos,allFragLengths,allDetectedMol,allMoleculesIm,timeInstances] = diffusion_process_image_speed(...
    cutInfo,constVal);



% Get the observed and true cumulative max time series 

timeSeriesObs = zeros(1,length(allDetectedMol));
nrVisCuts = timeSeriesObs;
prevNr = 0;


for i = 1:length(timeSeriesObs)	
	
	currNr = length(allDetectedMol{i});
	nrVisCuts(i) = currNr;
    
    if currNr > prevNr
        timeSeriesObs(i) = currNr;
        prevNr = currNr;
        
    else
        timeSeriesObs(i) = prevNr;
        
    end


end

timeSeriesTrue = zeros(1,length(allDetectedMol));
allCutTimes = cutInfo.cutTimes;

for i=1:length(timeInstances)
       
    currTime = timeInstances(i);
    
    % Loop through all cut times to find out which one we are located at
    for j=1:length(allCutTimes)
        
        timeSeriesTrue(i) = j;
        currCutTime = allCutTimes(j);
        
       if currCutTime > currTime
           break;
       end
        
    end
    
end



end

