function [cuttingInformation,timeToN0] = nicking_process_2(L,frayDist,numberOfCuts,alpha,initialNickNr)

% Needed functions
import nicking_simulation.check_cut;

% Create the original dna array
dna = zeros(L,2);

nrNonNicked = L*2;

% Create cell array to store all fragments in
fragments = cell(1);


% Insert the only fragment at the moment
fragments{1} = dna;


% Create cell array to store all the non nicked sites in
nonNickedSites = cell(1);
% Insert all the current non nicked sites from the dna array
[nonNickedSites{1}(:,1),nonNickedSites{1}(:,2)] = find(dna==0);

% Array to store cut times in
cutTimes = zeros(1,numberOfCuts);
cutAtFragNr = zeros(1,numberOfCuts);
cutPositions = zeros(1,numberOfCuts);
fragmentAfterCut = cell(1,numberOfCuts);
nrNicksAtCut = zeros(1,numberOfCuts);

% Array to store the curent number of fragments during the simulation in
% numberOfFrag = zeros(1,maxNickNr);
% numberOfSitesLeft = zeros(1,maxNickNr);

% Create time variable to know when cut happened
time = 0;
cutCount = 0;
timeToN0 = 0;
nickNr = 0;

% Loop through all the nickings
while cutCount < numberOfCuts
    
    nickNr = nickNr + 1;
    
    % Update the current time according to Gillespie
    time = time - log(rand)/(sum(nrNonNicked)*alpha);
    
    if nickNr == initialNickNr
        timeToN0 = time;
    end
    
    
    % Choose a random site to nick on a random fragment
    
    % Use the weighted random generator (bebbington) for choosing the fragment
    fragToNick = randsample(length(fragments),1,'true',nrNonNicked);
    
    % Non weigthed random generator for the nick site
    randPos = randi(nrNonNicked(fragToNick));
    
    % Extract the randomly choosen row and column on the fragment
    siteToNickRow = nonNickedSites{fragToNick}(randPos,1);
    siteToNickCol = nonNickedSites{fragToNick}(randPos,2);
    
    
    
    
    % Per construction we now know that the outcome is either a new nick or
    % a cut, so we first need to check for a cut before placing a nick
    
    % Check if cut happend
    [foundCut,cutPosRow] = check_cut(fragments{fragToNick},frayDist,siteToNickRow,siteToNickCol);
    
    if foundCut
        
        % Now we need to split the current fragment into two new ones and
        % conserve the spatial order of the fragments
        
        
        % Create the two new splited up fragments
        if max(siteToNickRow,cutPosRow) == min(siteToNickRow,cutPosRow)
            if max(siteToNickRow,cutPosRow) == size(fragments{fragToNick},1)
                secondFrag = fragments{fragToNick}(max(siteToNickRow,cutPosRow):end,:);
                firstFrag = fragments{fragToNick}(1:min(siteToNickRow,cutPosRow),:);
            else
                secondFrag = fragments{fragToNick}(max(siteToNickRow,cutPosRow)+1:end,:);
                firstFrag = fragments{fragToNick}(1:min(siteToNickRow,cutPosRow),:);
            end
        else
            secondFrag = fragments{fragToNick}(max(siteToNickRow,cutPosRow):end,:);
            firstFrag = fragments{fragToNick}(1:min(siteToNickRow,cutPosRow),:);
        end
        
        
        
        
        % New fragment cell array with one more element
        fragSubs = cell(1,length(fragments)+1);
        
        % Transfer all fragments to the new cell array preserving the
        % spatial order
        counter = 1;
        for fragNr=1:length(fragments)
            
            if fragNr == fragToNick
                fragSubs{counter} = firstFrag;
                counter = counter + 1;
                fragSubs{counter} = secondFrag;
                counter = counter + 1;
            else
                fragSubs{counter} = fragments{fragNr};
                counter = counter + 1;
            end
            
        end
        
        % Update the fragment list accordingly
        fragments = fragSubs;
        
        
        % Since cut happened we need to update the fragements
        
        % Array to store number of non-nicked sites in for all fragments
        nrNonNicked = zeros(1,length(fragments));
        
        % Find all the non nicked sites on all DNA fragments
        for fragNr=1:length(fragments)
            [rows,cols] = find(fragments{fragNr} == 0);
            nonNickedSites{fragNr} = [rows,cols];
            nrNonNicked(fragNr) = length(nonNickedSites{fragNr}(:,1));
        end
        
        
        % Increment the cut count varaible
        cutCount = cutCount + 1;
        % Save the position of the cut
        cutPositions(cutCount) = cutPosRow;
        % Save the fragment number at which the cut happened
        cutAtFragNr(cutCount) = fragToNick;
        % Save the time at which the cut happened
        cutTimes(1,cutCount) = time;
        % Save the number of nicks we accumulated at the current cut
        nrNicksAtCut(1,cutCount) = nickNr;
        % Save the current fragement costalation after the cut happened
        fragmentAfterCut{cutCount} = fragments;
        
        
        
    else
        % Now we just need to place a nick on the right site in the rigth
        % fragment
        fragments{fragToNick}(siteToNickRow,siteToNickCol) = 1;
        
        % And update the non nicked sites in the current fragment
        
        [rows,cols] = find(fragments{fragToNick} == 0);
        nonNickedSites{fragToNick} = [rows,cols];
        
        nrNonNicked(fragToNick) = nrNonNicked(fragToNick) - 1;
        
        
    end
    
    % Store the current number of fragments and total number of nickable
    % sites left
    %numberOfFrag(nickNr) = length(fragments);
    %numberOfSitesLeft(nickNr) = sum(nrNonNicked);
    
    
    if sum(nrNonNicked) == 0
        break % leave nicking process if no more positions can be nicked
    end
    
    
end

cuttingInformation.cutTimes = cutTimes;
cuttingInformation.cutPositions = cutPositions;
cuttingInformation.cuttedFragmentNr = cutAtFragNr;
cuttingInformation.NrNicksAtCut = nrNicksAtCut;
cuttingInformation.fragmentCollection = fragmentAfterCut;

end
