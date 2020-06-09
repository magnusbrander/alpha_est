%function produce_rate_equations(alpha,startInd,stopInd)

alpha = 2e-3
startInd = 1
stopInd = 4

% Specify needed constants 
constVal.L = 48490; %[basepairs]
constVal.pixelSize = 160; %[nanometers]
constVal.frayingDist = 1; %[basepairs]
constVal.diffConst = 628e3; %[nanometers^2 /s]
constVal.salt = 0.5;
constVal.intialNickDensity = 0;
constVal.timeResolution = 0.1;
constVal.waveLength = 509;
constVal.NA = 1.4;
constVal.uniSegConst = 5;
constVal.lambdaBg = 100;
constVal.lambdaSig = 1300;
constVal.gain = 254;


% Set expexted number of cuts to be simulated
expNrCuts = 40;

% Calculate the needed simulation time 
constVal.maxTime = sqrt(expNrCuts/(alpha^2*constVal.L * constVal.frayingDist));

% Set the number of time series to be simulated in order to create the empirical r(t) function 
statNr = abs(stopInd-startInd)+1;

import rate_function.empirical_rate_equation_image;


% Simulate the time series
simRateEq = empirical_rate_equation_image(alpha,constVal,statNr);

% Create filename 
filename = strcat('simulated_rate_equation_alpha_image',num2str(alpha),...
    '_',num2str(startInd),'_',num2str(stopInd),'.mat');

% Save file
save(['rate_equations/',filename],'simRateEq');


fprintf(strcat('Completed simulation for alpha=',num2str(alpha)))
fprintf('\n')



%end
