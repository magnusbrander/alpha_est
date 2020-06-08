function finalNumber = numInInterval(number,minNum,maxNum)

if number>minNum && number<maxNum
    finalNumber = number;
elseif number<minNum
    finalNumber = minNum;
else
    finalNumber = maxNum;
end

end

