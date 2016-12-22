function [goodsegs] = linesNearMask(I, M, lvl)
% given an image I and mask M, use Hough transformation to detect the lines
% that is near the mask rectangle.
% adjustable parameters: 
% threshold in houghpeaks
% fillgap and minlength in houghlines
% 

if ~ismatrix(I)
    I = rgb2gray(I);
end

% Gaussian filter before Hough
gaus = fspecial('gaussian', [5 5]);
I = imfilter(I, gaus, 'symmetric');

[R, C] = size(M);
[row, col] = find(M);
MRect = getRect(row, col); % four vertices of the rect
BW = edge(I,'log');
BW(M) = 0;
% BW = bwareaopen(BW, 2,4); % reduce small points
[H,T,R] = hough(BW, 'Theta', -90:0.5:89.5); % , 'RhoResolution', 0.5);
P = houghpeaks(H,6,'threshold',ceil(0.4*max(H(:)))); % the threshold should be reconsidered
fillgap = lvl*4;
expanddist = fillgap;
eSegThres = lvl*4;
minlength = lvl*3;
segs = houghlines(BW,T,R,P,'FillGap',fillgap,'MinLength',minlength); % FillGap, MinLength
nsegs = length(segs);
linest = [1];
goodsegs = [];

k = 1;
segs(1).k = 1;
for i = 2:nsegs
    if segs(i).rho ~= segs(i-1).rho || segs(i).theta ~= segs(i-1).theta
        linest = [linest, i];
        k = k+1;
    end
    segs(i).k = k;
end
linend = [linest(2:end)-1, nsegs];

% determine segments that is near the rectangle of mask
% isgoodseg = zeros(1, nsegs);
for i=1:nsegs
    segs(i).isgood = isGoodSeg(segs(i), MRect, expanddist);
end
i=1;
while(i<=nsegs)
    if segs(i).isgood>0
        s = linest(segs(i).k);
        e = linend(segs(i).k);
        goodsegs = [goodsegs, expandSeg(segs(s:e), i-s+1, eSegThres)];
        i = e;
    end
    i = i+1;
end

showHoughLines(I, M, goodsegs);
% close all;

end

function [eSeg] = expandSeg(segs, idx, thres)
% segs(idx) is a 'goodSeg' that have been found. Expand this segment so
% that it can be longer and 

pst = segs(idx).point1;
pend = segs(idx).point2;
len = length(segs);

% search pst
for i=idx-1:-1:1
    if segs(i).isgood || norm(pst - segs(i).point2) < thres
        pst = segs(i).point1;
    end
end

% search pend
for i=idx+1:len
    if segs(i).isgood || norm(pend - segs(i).point1) < thres
        pend = segs(i).point2;
    end
end

eSeg = segs(idx);
eSeg.point1 = pst;
eSeg.point2 = pend;
end

function [flag] = isGoodSeg(seg, rect, expandDist)
% if a segment cuts through the original rectangle and at least one vertice
% is in the expanded rectangle (from the original one by expandDist), it is
% considered a good Segment.

eRow = [rect.tl(1)-expandDist, rect.br(1)+expandDist];
eCol = [rect.tl(2)-expandDist, rect.br(2)+expandDist];
eRect = getRect(eRow, eCol);

p1 = seg.point1;
p2 = seg.point2;
p1 = [p1(2), p1(1)];
p2 = [p2(2), p2(1)];

flag = 0;
% diagonals are enough, four edges are redundant
if ~isCut(seg, rect.tl, rect.br) && ~isCut(seg, rect.bl, rect.tr) 
    return;
end
if ptInRect(p1, eRect) || ptInRect(p2, eRect) || diagCutSeg(seg, rect)
    flag = 1;
end

end

function [flag] = ptInRect(p, rect)

f1 = (p(1)-rect.tl(1)) * (rect.br(1)-p(1));
f2 = (p(2)-rect.tl(2)) * (rect.br(2)-p(2));
flag = f1>=0 && f2>=0;
end

function [flag] = diagCutSeg(seg, rect)
flag = 0;
p1 = seg.point1;
p2 = seg.point2;

% y = k1*x+b1 -> tl-br
x1 = rect.tl(2); y1 = rect.tl(1);
x2 = rect.br(2); y2 = rect.br(1);
k1 = (y2-y1)/(x2-x1);
b1 = y1 - k1*x1;

flag = (k1*p1(1)-p1(2)+b1)*(k1*p2(1)-p2(2)+b1)<=0;
if flag==1
    return;
end

% y = k2*x+b2 -> bl-tr
x1 = rect.bl(2); y1 = rect.bl(1);
x2 = rect.tr(2); y2 = rect.tr(1);
k2 = (y2-y1)/(x2-x1);
b2 = y1 - k2*x1;
flag = (k2*p1(1)-p1(2)+b2)*(k2*p2(1)-p2(2)+b2)<=0;
end

function [flag] = isCut(seg, A, B)

x1 = A(2);
y1 = A(1);
x2 = B(2);
y2 = B(1);
rho = seg.rho;
theta = seg.theta;
a = cosd(theta); b = sind(theta);

flag = (a*x1+ b*y1 -rho) * (a*x2+b*y2-rho)<=0;

end

function [MRect] = getRect(row, col)
rmin = min(row);
rmax = max(row);
cmin = min(col);
cmax = max(col);
MRect.tl = [rmin, cmin];
MRect.tr = [rmin, cmax];
MRect.bl = [rmax, cmin];
MRect.br = [rmax, cmax];

end