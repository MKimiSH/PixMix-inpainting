function [refI] = lightnessAdjust(refI, M, cont, lastI, lastcont)
% Lightness adjustment. 
% _cont_ and _lastcont_ should be corresponding contour points.
% Adjust the color value of the pixels of _refI_ in _M_ according to the
% difference of lightness between _cont_ and _lastcont_

[R, C, ~] = size(refI);
hst = 8; hszg = hst+1;
vst = 8; vszg = vst+1;
sizegrid = hszg*vszg;
% rI = refI;

%% 计算边缘点的像素值以及差值
% 我认为不需要用所有边缘点。
nborder = length(cont);
rind = sub2ind(cont(:,1), cont(:,2)); % which is x, which is y?
lind = sub2ind(lastcont(:,1), lastcont(:,2));
tl = min(cont);
br = max(cont);

b = zeros(nborder, 3);
bpr = zeros(nborder, 3);
b(:,1) = refI(rind); b(:,2) = refI(rind+R*C); b(:,3) = refI(rind+R*C*2);
bpr(:,1) = lastI(lind); bpr(:,2) = lastI(lind+R*C); bpr(:,3) = lastI(lind+R*C*2);
bdiff = b - bpr;

%% 计算网格G以及\mu(g_i)
xgv = round(linspace(tl(2), br(2), hst+1));
ygv = round(linspace(tl(1), br(1), vst+1));
[X, Y] = meshgrid(xgv, ygv);
colgrid = X(:); 
rowgrid = Y(:);
g = [rowgrid, colgrid];
dist = abs(bsxfun(@minus, rowgrid', cont(:,1)) + ...
           bsxfun(@minus, colgrid', cont(:,2)));

expdist = exp(-sqrt(dist));
thetagi = sum(expdist);
mugi1 = sum(bsxfun(@times, bdiff(:,1), expdist))./thetagi;
mugi2 = sum(bsxfun(@times, bdiff(:,2), expdist))./thetagi;
mugi3 = sum(bsxfun(@times, bdiff(:,3), expdist))./thetagi;

%% 利用插值计算矩形(tl, br)之中每个点p的\mu(p)值
% 这些3行的东西应该可以改成循环或者向量化，但是我懒得>_<
% [Xq, Yq] = meshgrid(tl(2):br(2), tl(1):br(1));
Xq = tl(2):br(2);
Yq = tl(1):br(1);
V1 = reshape(mugi1, [hszg, vszg]);
V2 = reshape(mugi2, [hszg, vszg]);
V3 = reshape(mugi3, [hszg, vszg]);
Vq1 = interp2(xgv, ygv, V1, Xq, Yq);
Vq2 = interp2(xgv, ygv, V2, Xq, Yq);
Vq3 = interp2(xgv, ygv, V3, Xq, Yq);
adjmat = zeros(R, C, 3);
adjmat(Yq, Xq, 1) = Vq1;
adjmat(Yq, Xq, 2) = Vq2;
adjmat(Yq, Xq, 3) = Vq3;
mind = find(M>0);
refI(mind) = refI(mind) + adjmat(mind);
refI(mind+R*C) = refI(mind+R*C) + adjmat(mind+R*C);
refI(mind+R*C*2) = refI(mind+R*C*2) + adjmat(mind+R*C*2);
end