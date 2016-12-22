function oddNodes = pointInPolygon(xy, array)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    l = size(array,1);
    j = l - 1;
    x = xy(1);
    y = xy(2);
    oddNodes = 0;
    for i = 1:l
        if((array(i,2) < y && array(j,2) >= y   ...
            || array(j, 2) < y && array(i, 2) >= y) ...
            && (array(i,1) <= x || array(j, 1) <= x))
            oddNodes = xor( oddNodes, (array(i,1)+(y - array(i,2))... 
                /(array(j,2)-array(i,2))*(array(j,1)-array(i,1))<x));
        end
        j = i;    
    end
end

