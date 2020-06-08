function nested_sampling()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference 
% for the data specified in
% misc.data_id together with data path. A number of parameters
% may be specified by the user:
%
%  - Dmin, Dmax, mumin, mumax
%    specify the minimum and maximum values of the two parameters   
%    which model the dynamics, and must be chosen to ensure that 
%    the most likely volume of parameter space is within these ranges.
%
%  - options.nwalkers is the number of initial samples used in the
%    nested sampling algorithm to sample parameter space.
%
%  - options.stoprat is the ratio between the last added evidence and the
%    total samples evidence at which nested sampling will terminate.
%
%  - options.nsteps is the number steps attempted with the MCMC
%    technique to generate a uniformly distributed sample from another sample.
%
%  - the "models" structure specify the functions used for likelihood calculation,
%    sample generation, which are specific 
%    to the system of interest. To use the nested sampling framework on another system, 
%    similar functions must exist and be specified in this structure.
%
% The routine outputs a .mat file with information on the user set parameters
% as well as evidence estimations and inferred parameter means and the samples used.
% In addition a .txt file is written, holding the conclusions of the analysis.
% The routine writing this file takes the following inputs:
%  - misc.percentiles are the percentiles used for the characterization of the
%    posterier 
%  - misc.labels are the labels assigned to each inferred parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


% Specify constants (guess)
L = 48000;
Xi = 2;

%Specify a uniform prior on alpha and n0
alphamax=10^(-1);   
models(1).invprior=@(u) [u(1)*alphamax; u(2)];

%Specify options
options.nwalkers=200;   % Number of walkers to be generated
options.stoprat=10^(-3);% Ratio to stop sampling
options.nsteps=30;      % Attempted number of steps in parameter space
models(1).options=options;

%Specify the u-generators
models(1).genu=@() rand(1,2);

%Specify the logl
import timeseries_folder.cutProb;
models(1).logl=@(obs,theta) -cutProb(obs.obs,theta(1),obs.L,Xi,theta(2));

%Specify the index for the labels
models(1).labels=[1 2];

%Percentiles to be calculated
misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];

%Labels for the parameters and percentile line in the output text-file
misc.labels=...
['alpha: ';...
 'n0:    '];

%Specify output filename beginnings
misc.data_id = 'simple';


%I = csvread('/home/magnus/Documents/MATLAB/SDD_modified/all_timeSeries_5x_oxygen.txt');
%framesPerSec = 10;
I = csvread('/home/magnus/Documents/MATLAB/SDD_modified/timeSeries_simulation.txt');
framesPerSec = 1;
nrCuts = 4;

ITrimed = zeros(1,nrCuts);

count = 1;
for i = 1:size(I,1)
    
    iTemp = I(i,:);
    iTemp = iTemp(iTemp ~= 0);
    if length(iTemp) >= nrCuts
        
        ITrimed(count,:) = iTemp(1:nrCuts);
        count = count+1;
    end
        
    
end


lengthMatrix = L * ones(size(ITrimed,1),size(ITrimed,2));

timeOfCuts = cell(1,size(ITrimed,1));
allLengths = timeOfCuts;


for i=1:size(ITrimed,1)
    
    timeOfCuts{i} = ITrimed(i,:)/(framesPerSec);
    allLengths{i} = lengthMatrix(i,:);
    
end

data.obs = timeOfCuts;
data.L = allLengths;


%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Run the nested sampling algorithm for all models and compile results
addpath('/home/magnus/Documents/MATLAB/SDD_modified/nested_sampling/nestedsampling-master');
[results] = ns_processdataset(data,models,misc);

path=[misc.data_id,'_output'];

%Save results
save(path,'results')
save(path,'data','-append')
