function [cutTimes,nrNicksAtCut,timeToN0] = nicking_process(L,frayDist,numberOfCuts,alpha,initialNickNr)

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
nrNicksAtCut = zeros(1,numberOfCuts);

% Array to store the curent number of fragments during the simulation in
%numberOfFrag = zeros(1,maxNickNr);
%numberOfSitesLeft = zeros(1,maxNickNr);

% Create time variable to know when cut happened
time = 0;
cutCount = 0;
timeToN0 = 0;
nickNr = 0;

% Loop through all the nickings 
while cutCount < numberOfCuts
    
    nickNr = nickNr + 1;
    
    % Update the current time at which the nick/cut happened
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
    % a cut so we first need to check for a cut before placing a nick   
    
    % Check if cut happend 
    [foundCut,cutPosRow] = check_cut(fragments{fragToNick},frayDist,siteToNickRow,siteToNickCol);
    
    if foundCut
        
        % Now we need to split the current fragment into two new ones
        
        % We can do this by adding the second part of the current fragment
        % as a new fragement and after that replace the current fragment
        % with the first part of the current fragment 
        
        secondFrag = fragments{fragToNick}(max(siteToNickRow,cutPosRow):end,:);
        firstFrag = fragments{fragToNick}(1:min(siteToNickRow,cutPosRow),:);
        
        fragments{end+1} = secondFrag;
        fragments{fragToNick} = firstFrag;
        
        
        % If cut happened we need to update the fragements
        
        % Array to store number of non-nicked sites in for all fragments
        nrNonNicked = zeros(1,length(fragments));
        
        % Find all the non nicked sites on all DNA fragments
        for frag=1:length(fragments)
            [rows,cols] = find(fragments{frag} == 0);
            nonNickedSites{frag} = [rows,cols];
            nrNonNicked(frag) = length(nonNickedSites{frag}(:,1));
        end
        
        
        cutCount = cutCount + 1;
        cutTimes(1,cutCount) = time;
        nrNicksAtCut(1,cutCount) = nickNr;
        
        
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
        return
    end
    
 
end




end
