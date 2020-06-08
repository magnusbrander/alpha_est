function currPositions = assign_new_pos(prevLengths,currLenghts,prevPositions,cuttedFragNr)
% FUNCTION: newPos = assign_new_positions(
% This function places the fragments in the right positions after a cut
% happened taking into account the previous positions of the fragments




        currPositions = zeros(1,length(currLenghts));
        counter = 1;
        for k=1:length(prevLengths) % go through all fragments in previous instance
            
            if k == cuttedFragNr
                
                % Assign positions to the new fragments based on the
                % position of the parent fragment
                firstNewPos = prevPositions(k) + 0.5*(currLenghts(counter) -prevLengths(k));
                currPositions(counter) = firstNewPos;
                counter = counter +1;
                secondNewPos = prevPositions(k) + 0.5*(prevLengths(k) - currLenghts(counter));
                currPositions(counter) = secondNewPos;
                counter = counter +1;
                
            else
                
                currPositions(counter) = prevPositions(k);
                counter = counter + 1;
                
            end


        end


end
