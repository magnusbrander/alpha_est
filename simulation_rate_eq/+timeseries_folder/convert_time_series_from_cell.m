function arrayTimeSeries = convert_time_series_from_cell(cellTimeSeries,minCutNr)

% ###################################################################
%
% FUNCTION: arrayTimeSeries = convert_time_series_from_cell(cellTimeSeries,minCutNr)
%
% Description: This converts a number of time series returned in a cell
% array format to an array list 
%
% Input: cellTimeSeries = all time series in cell array format, minCutNr =
% the number of cutting times that should be included
%
% Output: arrayTimeSeries = an array with all the indivdual time series 
% listed on a separate row
%
% ###################################################################




% First we go throug all time series to check how many fullfill the
% recuired length criterion 

count = 0;

for i=1:length(cellTimeSeries)
    
    if length(cellTimeSeries{i}) >= minCutNr
        count = count + 1;
    end
        
end



% Now we can create the array to store all time series in

arrayTimeSeries = zeros(count,minCutNr);

% Here we transfer all the time series of sufficient length to the array
count = 0;
for i=1:length(cellTimeSeries)
    
    if length(cellTimeSeries{i}) >= minCutNr
        count = count + 1;
        arrayTimeSeries(count,:) = cellTimeSeries{i}(1:minCutNr);
    end
    
end



end
