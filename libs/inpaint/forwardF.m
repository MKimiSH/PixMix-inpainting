function [refI, initF] = forwardF(I, M, lI, lF, lM, H)
% Forward the mapping _lF_ to current frame _I_ and create a reference
% model _refI_ with homography _H_ 
% I worry about the precision loss caused by ROUND()...

[R, C] = size(lM);
assert(R == size(M,1) && C == size(M,2));

[crow, ccol] = find(M>0);
cn = length(crow);
%% reverse transform to lI
currP = [ccol, crow, ones(cn,1)]; % row->y, col->x
refP = currP*(inv(H)); % row vectors
rrow = round(refP(:,2)./refP(:,3));
rcol = round(refP(:,1)./refP(:,3));
%% find mapping in lI
lind = sub2ind([R,C], rrow, rcol);
lfr = lF(lind);
lfc = lF(lind+R*C);
%% forward transform to I
refP = single([lfc, lfr, ones(cn,1)]);
currP = refP*H;
cfr = round(currP(:,2)./currP(:,3));
cfc = round(currP(:,1)./currP(:,3));
%% fill initF with mapping values
initF = zeros(R,C,2, 'int32');
[initF(:,:,2), initF(:,:,1)] = meshgrid(1:C, 1:R);
cind = sub2ind([R,C], crow, ccol);
initF(cind) = cfr;
initF(cind + R*C) = cfc;

%% Create refI
refI = I;
tform = projective2d(H);
warpLI = imwarp(lI, tform);
[WR, WC, ~] = size(warpLI);
t1 = [1,1,1]*H;
t1 = round(t1/t1(3));
lastrow = crow + t1(2); lastcol = ccol - t1(1);
lastind = sub2ind(size(warpLI), lastrow, lastcol);
refI(cind) = warpLI(lastind);
refI(cind + R*C) = warpLI(lastind + WR*WC);
refI(cind + 2*R*C) = warpLI(lastind + 2*WR*WC);
% figure, imshow(refI);

end