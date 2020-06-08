function numOut = appendNum(num1,num2)


% This function appends the num2 to num1 by the logic of: appendNum(3,1) ->
% 31. The output is a number and required input has to be in integer form. 

numOut = str2num(strcat(num2str(num1),num2str(num2)));

end
