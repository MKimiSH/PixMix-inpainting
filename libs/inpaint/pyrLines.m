function pyLines = pyrLines(pyI, pyM, L)

pyLines = cell(L, 1);

proclvl = L-1;
pyLines{proclvl} = linesNearMask(pyI{proclvl}, pyM{proclvl}, proclvl);
% showHoughLines(pyI{proclvl}, pyM{proclvl}, pyLines{proclvl});
nlines = length(pyLines{proclvl});
toplines = pyLines{proclvl};

for l = L:-1:1
    if l == proclvl
        continue;
    end
    fac = 2^(proclvl-l);
    for i = 1:nlines
        curline.point1 = (toplines(i).point1 - 1)/fac + 1; % original is (1,1)
        curline.point2 = (toplines(i).point2 - 1)/fac + 1;
        curline.theta = toplines(i).theta;
        curline.rho = toplines(i).rho/fac;
        pyLines{l} = [pyLines{l}, curline];
    end
end


end