function output = parfor_test_func(input)



% Take the input value and compute the squareroot

import timeseries_folder.centOfMass
output = centOfMass([input,input*2],[input,input*2]);



end