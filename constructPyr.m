function [pyI, pyM] = constructPyr(I, M)
% the rule is that in the first level, any px in the hole should be at
% most 3 px in distance from the 

% M should be logical to prevent interpolation in imresize()!!

Mcur = M;
L = 0;
for i=1:10 % 10 should be enough
    D = bwdist(~Mcur);
    maxDist = max(D(:));
    if maxDist <= 3
        L = i;
        break;
    end
    Mcur = imresize(Mcur, 0.5);
end

pyI = cell(L, 1);
pyM = cell(L, 1);

pyI{L} = I;
pyM{L} = M;
for l=L-1:-1:1
    pyM{l} = imresize(pyM{l+1}, 0.5);
    pyI{l} = imresize(pyI{l+1}, [size(pyM{l},1), size(pyM{l},2)] );
end

end