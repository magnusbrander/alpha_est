function [f,fittedParm] = fit_erf(curve,varargin)



nrOfArg = length(varargin);


 A0 = min(curve);
 F10 = max(curve)-min(curve);
 
 

 xData = 1:length(curve);
 yData = curve;
 

switch nrOfArg
    
    case 3
        
        % For one molecule
        xMax = length(curve);
        parmSumStart = 4;
        % Get input guesses
        [S0,B0,C0] = varargin{:};
        
        % Set starting points
        p0 = [A0, F10, S0, B0, C0];
        
        % Set lower and upper bounds 
        ALower = 0; AUpper = max(yData);
        BLower = 0; BUpper = xMax-1;
        CLower = 1; CUpper = xMax;
        SLower = S0; SUpper = S0;
        F1Lower = 0; F1Upper = Inf; 
        
        LB = [ALower, F1Lower, SLower, BLower, CLower];
        UB = [AUpper, F1Upper, SUpper, BUpper, CUpper];
        
        % Define the function to fit 
        % Define function to fit 
        y = @(p) p(1) + p(2) * ( erf(  (xData-p(4)) / ( sqrt(2)*p(3)) )  -...
            erf( (xData-(p(4)+p(5))) / (sqrt(2)*p(3))  )   );
        
        
        
    case 7
        
        % For two molecules
        
        % Get input guesses
        [S0,B0,C0,B20,C20,xMax,xMin] = varargin{:};
        
        parmSumStart = 5;
        
        F20 = F10;
        % Starting point for the model fitting.
        p0 = [A0 , F10, F20, S0, B0, B20, C0, C20];
        
        % Lower and upper constraints for the fitting
        ALower = 0; AUpper = max(yData);
        BLower = xMin; BUpper = xMax;
        B2Lower = 2; B2Upper = xMax;
        CLower = 2; CUpper = xMax;
        C2Lower = 2; C2Upper = xMax;
        SLower = S0; SUpper = S0;
        F1Lower = 0; F1Upper = Inf;
        F2Lower = 0; F2Upper = Inf;
        
        LB = [ALower , F1Lower, F2Lower, SLower, BLower ,B2Lower, CLower, C2Lower ];
        UB = [AUpper , F1Upper, F2Upper, SUpper, BUpper ,B2Upper, CUpper , C2Upper];
        
        % Define the function to fit
        y = @(p) p(1) + p(2) * (erf(  (xData-p(5)) / ( sqrt(2)*p(4)) )  -...
            erf( (xData-(p(5)+p(6))) / (sqrt(2)*p(4))  ) ) + ...
            p(3) * ( erf(  (xData- (p(5)+p(6)+p(7)) ) / ( sqrt(2)*p(4)) ) - ...
            erf( (xData-(p(5)+p(6)+p(7)+p(8))) / (sqrt(2)*p(4))  ) ) ;
        
        
    case 9
        
        % For 3 molecules
        
        % Get input guesses
        [S0,B0,C0,B20,C20,B30,C30,xMax,xMin] = varargin{:};
        
        parmSumStart = 6;
        
        F20 = F10;
        F30 = F10;
        % Starting point for the model fitting.
        p0 = [A0, F10, F20, F30, S0, B0, B20, B30, C0, C20, C30];
        
        % Lower and upper constraints for the fitting
        ALower = 0; AUpper = max(yData);
        BLower = xMin; BUpper = xMax;
        B2Lower = 1; B2Upper = xMax;
        B3Lower = 1; B3Upper = xMax;
        CLower = 1; CUpper = xMax;
        C2Lower = 1; C2Upper = xMax;
        C3Lower = 1; C3Upper = xMax;
        SLower = S0; SUpper = S0;
        F1Lower = 0; F1Upper = Inf;
        F2Lower = 0; F2Upper = Inf;
        F3Lower = 0; F3Upper = Inf;
        
        LB = [ALower ,F1Lower,F2Lower,F3Lower, SLower, BLower ,B2Lower, B3Lower, CLower, C2Lower, C3Lower ];
        UB = [AUpper ,F1Upper,F2Upper,F3Upper, SUpper, BUpper ,B2Upper, B3Upper, CUpper , C2Upper, C3Upper];
        
        
        % Define the function to fit
        y = @(p) p(1) + p(2) * ( erf(  (xData-p(6)) / ( sqrt(2)*p(5)) )  -...
            erf( (xData-(p(6)+p(7))) / (sqrt(2)*p(5))  ) )  + ...
            p(3) * ( erf(  (xData- (p(6) + p(7) + p(8)) ) / ( sqrt(2)*p(5)) ) - ...
            erf( (xData-(p(6)+p(7)+p(8)+p(9))) / (sqrt(2)*p(5))  ) ) + ...
            p(4) * (erf(  (xData- (p(6)+p(7)+p(8)+p(9)+p(10)) ) / ( sqrt(2)*p(5)) ) - ...
            erf( (xData-(p(6)+p(7)+p(8)+p(9)+p(10)+p(11))) / (sqrt(2)*p(5))  ) );
        
end



options = optimset('MaxIter',5000);

% Define objective as chi squared
objective = @(p) sum( (yData - y(p) ).^2 );


% Minimize the object with respect to the parameters
fittedParm = fmincon(objective,p0,[],[],[],[],LB,UB,@(p)constraints(p,xMax,parmSumStart),options);


% Get function
f = getFunc(fittedParm);



end






function [c,ceq] = constraints(p,xMax,parmSumStart)


% We require the sum of postions of the erf to be smaller than a certain
% value 
c = sum(p(parmSumStart:end)) - xMax;

% No inequality constraint is used
ceq = [];


end


function f = getFunc(parameters)

p = parameters;

switch length(p) 
    
    case 5
                
        f = @(x) p(1) + p(2) * ( erf(  (x-p(4)) / ( sqrt(2)*p(3)) )  -...
            erf( (x-(p(4)+p(5))) / (sqrt(2)*p(3))  )   );
        
    case 8
        
        f = @(x) p(1) + p(2) * (erf(  (x-p(5)) / ( sqrt(2)*p(4)) )  -...
            erf( (x-(p(5)+p(6))) / (sqrt(2)*p(4))  ) ) + ...
            p(3) * ( erf(  (x- (p(5)+p(6)+p(7)) ) / ( sqrt(2)*p(4)) ) - ...
            erf( (x-(p(5)+p(6)+p(7)+p(8))) / (sqrt(2)*p(4))  ) ) ;
                
    case 11
        
        f = @(x) p(1) + p(2) * ( erf(  (x-p(6)) / ( sqrt(2)*p(5)) )  -...
            erf( (x-(p(6)+p(7))) / (sqrt(2)*p(5))  ) )  + ...
            p(3) * ( erf(  (x- (p(6) + p(7) + p(8)) ) / ( sqrt(2)*p(5)) ) - ...
            erf( (x-(p(6)+p(7)+p(8)+p(9))) / (sqrt(2)*p(5))  ) ) + ...
            p(4) * (erf(  (x - (p(6)+p(7)+p(8)+p(9)+p(10)) ) / ( sqrt(2)*p(5)) ) - ...
            erf( (x -(p(6)+p(7)+p(8)+p(9)+p(10)+p(11))) / (sqrt(2)*p(5))  ) );
        
end


end



