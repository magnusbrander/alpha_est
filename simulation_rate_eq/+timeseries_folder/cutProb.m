function P = cutProb(timeOfCuts,alpha,L,Xi,n0)


% This function computes the probability for a series of
% time instances given a set of parameters 


% The "timeOfCuts" variable should be given as a cell array 
% The remaining arguments can be given as a number or a vector, but only
% one at a time can be a vector 





 % The likelihood function is given by:
% P(t_i|theta_j) = r(t_1)*r(t_2),...r(t_m) * exp(- (2*alpha*L*Xi*[alpha*T^2+T*n0]))


% r(t) is the rate of cuts at time t + dt , alpha the nicking rate,
% L the length of the fragment, n0 the number of initial nicks
% and T the total time 


% Extract the time of last cut and use this as the last trusted time
% instance

T = zeros(1,length(timeOfCuts));

for i=1:length(timeOfCuts)
    
    T(i) = timeOfCuts{i}(end);
    
end



% Create the likelihood log(P(t|theta)) as zero to begin with
P = 0;


% Compute the sum of the log likelihood functions including all time series
% and add them togheter

for seriesNr=1:length(timeOfCuts)
    
    for cutNr=1:length(timeOfCuts{seriesNr})
         
         
        P = P + log( 2*L{seriesNr}(cutNr)*alpha.* (1-exp(-(Xi)*(alpha*timeOfCuts{seriesNr}(cutNr) + n0)) ) );
        
        
    end
    
%     Add the balancing factor ln(exp(-
%     (2*alpha*L*Xi*[0.5*alpha*T^2+T*n0]))) to obtain the final likelihood for the given time series
     
    P =  P - ( 2 * mean(L{seriesNr}) * ( alpha*T(seriesNr) + (exp(-(Xi)*n0)./(Xi)).*(exp(-(Xi)*alpha*T(seriesNr) ) - 1 ) ) );

end


% Return the negative log likelihood
P = -P;



end
