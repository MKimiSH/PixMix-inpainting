function rgb = mf_smallst(m1, m2, i, j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    rgb = reshape((m1(i,j,:)*2 + m2(i,j,:))/3,1,3);

end

