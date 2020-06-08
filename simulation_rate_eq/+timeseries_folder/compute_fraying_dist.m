function [frayDist,constants] = compute_fraying_dist(T,C)
% #######################################################################
% This function computes the theoretical fraying distance assuming random
% and equal amounts of AT to GC 
%
% Input: 
% T = the temperature in celsius
% C = concentration of NaCl 
%
% Output: 
% fraying distance = twice the largest distance between nicks to get a cut 
% #######################################################################



%Input temperature in celsius 

T = 273.15 + T;

Tref = 273.15+37;

% Compute the RT factor 
RT = 1.9859*10^(-3)*T;


% Salt concentration 
Cref = 0.1;

% Compute melting temperature for the different bases depndent on NaCl
% concentration
Tm_AT = 355.55 + 7.95*log(C);
Tm_GC = 391.55 + 4.89*log(C);


% Values for geometric free energy cost at 37 C and 0.1M NaCl
GrefATTA = - 1.73; 
GrefTAAT = - 0.58;  
GrefAATT = - 1.50;
GrefGACT = - 1.82; 
GrefCAGT = - 0.94; 
GrefAGTC = - 1.45; 
GrefACTG = - 2.20;
GrefGGCC = - 1.83;
GrefCGGC = - 1.30;
GrefGCCG = - 2.56;


% Compute stacking energies corrected for temperature and salt
% concentration 

% Correction for salt 
GrefATTA_salt = GrefATTA - 0.205*log(C/Cref);
GrefTAAT_salt = GrefTAAT - 0.205*log(C/Cref);
GrefAATT_salt = GrefAATT - 0.205*log(C/Cref);
GrefGACT_salt = GrefGACT - 0.205*log(C/Cref);
GrefCAGT_salt = GrefCAGT - 0.205*log(C/Cref);
GrefAGTC_salt = GrefAGTC - 0.205*log(C/Cref);
GrefACTG_salt = GrefACTG - 0.205*log(C/Cref);
GrefGGCC_salt = GrefGGCC - 0.205*log(C/Cref);
GrefCGGC_salt = GrefCGGC - 0.205*log(C/Cref);
GrefGCCG_salt = GrefGCCG - 0.205*log(C/Cref);

% Correction for temperature
G_st_temp = 0.026*(T-Tref);

% Corrected stacking energies thus become
Gst_ATTA = GrefATTA_salt + G_st_temp;
Gst_TAAT = GrefTAAT_salt + G_st_temp;  
Gst_AATT = GrefAATT_salt + G_st_temp;
Gst_GACT = GrefGACT_salt + G_st_temp; 
Gst_CAGT = GrefCAGT_salt + G_st_temp; 
Gst_AGTC = GrefAGTC_salt + G_st_temp; 
Gst_ACTG = GrefACTG_salt + G_st_temp;
Gst_GGCC = GrefGGCC_salt + G_st_temp;
Gst_CGGC = GrefCGGC_salt + G_st_temp;
Gst_GCCG = GrefGCCG_salt + G_st_temp;




% Free energy associated with hydrogen bond for AT and GC

% Two cases for AT and GC repectively assuming 50:50 AT:GC ratio


% For the AT case
deltaG_hb_AT = -24.85*10^(-3)*(Tm_AT - T) - 0.25*(Gst_ATTA+Gst_TAAT+2*Gst_AATT);

% For the GC case
deltaG_hb_GC = -24.85*10^(-3)*(Tm_GC- T) - 0.25*(Gst_GCCG+Gst_CGGC+2*Gst_GGCC);

deltaG_hb = 0.5*(deltaG_hb_AT + deltaG_hb_GC);




% Here we do not take any informatio about the neighboring basepairs so we
% compute an average of the stacking energies to resemble the final
% contribution 

avgGst =  (Gst_ATTA + Gst_TAAT + 2*Gst_AATT + 2*Gst_GACT + 2*Gst_CAGT + 2*Gst_AGTC +...
    2*Gst_ACTG + 2*Gst_GGCC + Gst_CGGC + Gst_GCCG)/16;

% The final energy coast for opening a base pair is then given by 
deltaG = deltaG_hb + 2*avgGst;
 
r = exp( (avgGst + deltaG_hb)/RT);
k = exp(avgGst/RT); 
% With the free energy cost we estimate the fraying distance as 
%frayDist = 2 * (exp(deltaG./RT)./(1-exp(deltaG./RT))) + 1;
frayDist = 2 * (r*k)/((1-r)*(1-r+k*r));
% The extra one comes from the fact that we have cut for all values below 1
% due to the basepair configuration

constants.Gst = avgGst;
constants.Ghb = deltaG_hb;
constants.RT = RT;


end
