function logP = diffusion_log_likelihood(allPosSeries,D,dt)



% #########################################################################
%
% FUNCTION: output = diffusion_log_likelihood(input)
%
%
% Description: 
%
% Input: 
%
% Output: 
%
% #########################################################################






logP = 0;

constLog = log(sqrt(4*pi*D*dt));


for i = 1:length(allPosSeries)
    
    x = allPosSeries{i};
    T = length(allPosSeries{i});
    
    for t = 2:1:T
        
        logP = logP - ((x(t)-x(t-1))^2)/(4*D*dt) - constLog;
        
    end
    
    
    
end

end
