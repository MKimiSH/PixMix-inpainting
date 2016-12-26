function pyLines = pyrLines(pyI, pyM, usrLines, L)

pyLines = cell(L, 1);

proclvl = L-1;
pyLines{proclvl} = linesNearMask(pyI{proclvl}, pyM{proclvl}, proclvl);
% showHoughLines(pyI{proclvl}, pyM{proclvl}, pyLines{proclvl});
nlines = length(pyLines{proclvl});
toplines = pyLines{proclvl};

%% 去个重
if ~isempty(usrLines)
    fac = (proclvl-L);
    nusrln = length(usrLines);
    fprintf('%d lines detected, user added %d lines\n', nlines, nusrln);
    for i = 1:nusrln
        usrLines(i).point1 = (usrLines(i).point1 - 1)/fac + 1; % original is (1,1)
        usrLines(i).point2 = (usrLines(i).point2 - 1)/fac + 1;
        usrLines(i).theta = usrLines(i).theta;
        usrLines(i).rho = usrLines(i).rho/fac;
        for j = 1:nlines
            if linesTooClose(usrLines(i), toplines(j))
%                 toplines(j) = usrLines(i); 
                % too close: trust the detected line
                break;
            end
        end
        if isempty(j) || j == nlines % no close: add the user defined line
            toplines = [toplines, usrLines(i)];
        end
    end
end

%% 扩到每一个level
nlines = length(toplines);
fprintf('%d lines total.\n', nlines);
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

function [isclose] = linesTooClose(l1, l2)
if l1.theta <0
    l1.theta = l1.theta+180;
    l1.rho = -l1.rho;
end
if l2.theta <0
    l2.theta = l2.theta+180;
    l2.rho = -l2.rho;
end

isclose = (abs(l1.theta - l2.theta)<2) && (abs(l1.rho-l2.rho)<2);

end