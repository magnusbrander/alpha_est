function [ratesLeft,ratesRight] = get_rates(fragmentLengths,D_lambda,L,p)

% FUNCTION: rates = get_rates(fragmentLengths)
% This functions intends to compute the hop rate for each fragment based on
% their individual length assuming the diffusion constant scales as 1/L
% Where L is the orignal length




% Convert rates for the  other diffusion constant based on the different lengths 
 ratesLeft = (D_lambda/p^2)*(L./fragmentLengths);
 ratesRight = ratesLeft;




end
