function S = poissrnd_speed(lambda,dim1,dim2)

% Generate a random sample S of size ns from the (discrete)
% Poisson distribution with parameter lambda.
% Fixed error:
%    CHANGED k = 1; produ = 1; produ = produ*rand
%    TO      k = 0; produ = rand;
% Derek O'Connor, 24 July 2012.  derekroconnor@eircom.net
%
S = zeros(dim1,dim2);
for row = 1:dim1
    
    for col=1:dim2
        
        k = 0;
        produ = rand;
        while produ >= exp(-lambda)
            produ = produ*rand;
            k = k+1;
        end
        S(row,col) = k;
        
    end
    
end