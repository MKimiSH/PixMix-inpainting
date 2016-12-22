function rgb = mf(m1, m2, i, j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    rgb = reshape(uint8((double(m1(i,j,:))*2 + double(m2(i,j,:)))/3),1,3);

end

