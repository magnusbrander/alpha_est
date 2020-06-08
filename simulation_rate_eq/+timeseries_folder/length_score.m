function lenMatchIndex = length_score(prevLen,currLengths,indOfMatch)


% Magic number based on variance

const = 0.2;

tollerance = 1.5 * const*prevLen;

lenDiff = abs(currLengths-prevLen);

lenMatchIndex = zeros(1,length(indOfMatch));

for i=1:length(indOfMatch)
    
    if lenDiff(indOfMatch(i))<tollerance
        lenMatchIndex(i) = 1;
    end
end

end


