function gradVec = modified_gradient(y,dx)


% ########################################################################
%
% FUNCTION: gradVec = modified_gradient(y,dx)
%
% Description: Takes the central difference of a vector while truncating
% the difference in the ends to only forward and backward in the begining
% and end, respectively.
%
% Input: y = vector of values, dx = step size between values
%
% Output: gradVec = central difference of y 
%
% ########################################################################


gradVec = 0*y;


for i=1:length(y)
    
    if i==1
        
        gradVec(i) = (y(i+1) - y(i))/dx;
        
    elseif i==length(y)
        
        gradVec(i) = (y(i) - y(i-1))/dx;
        
    else
        
        gradVec(i) = (y(i+1) - y(i-1))/(2*dx); 
        
    end
    
    
end


end

