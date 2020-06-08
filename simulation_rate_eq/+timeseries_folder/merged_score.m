function [intactScore,matchIndex] = merged_score(totFragToTrack,potMatchFrag,image,time)


% This function intend to compare to fragments and see if they have merged
% into another one 

% We first want to create ONE new fragment from the two input fragments in
% order to use the intact_score function 

% Estimate the length of the fragments
firstFragLenEst = sum(totFragToTrack(:,2))/length(totFragToTrack(:,2));
secondFragLenEst = sum(potMatchFrag(:,2))/length(potMatchFrag(:,2));
lengths = [firstFragLenEst,secondFragLenEst];


% The total length of the new fragment is estimated to be the sum of their
% individual estimated lengths
totLengthEst = firstFragLenEst + secondFragLenEst;

% Here we compute their joint center of mass
prevPosFirst = totFragToTrack(size(totFragToTrack,1),1);
prevPosSecond = potMatchFrag(size(potMatchFrag,1),1);

positions = [prevPosFirst,prevPosSecond];

% Here we compute the two fragments weigthed center of mass position. This
% value does not need to be rounded
import timeseries_folder.centOfMass
jointCentMass = centOfMass(positions,lengths);


% Create the substitute fragment 
subFragment = [jointCentMass,totLengthEst,0,0];

% Here we use the intact_score function to check if our hypotesis that
% these two fragments merged is consistent with the current fragments 

import timeseries_folder.intact_score
% Here we return a logical value answering yes or no depending on if it is intact 
% as well as which fragment index that gave the potential match in the
% current time instance
[intactScore,matchIndex] = intact_score(subFragment,image,time);



end