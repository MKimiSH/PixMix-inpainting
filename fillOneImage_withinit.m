function [outI, outF] = fillOneImage_withinit(I, M, prevI, prevMap)
if ~ismatrix(M)
    M = rgb2gray(M);
end
I = im2single(I);
M = im2single(M);
M = M>0;

I(M==1) = 0;

tst = tic;

figure; imshow(I);
figure; imshow(prevI);

prevMap = int32(prevMap);
H = estimateTransform(I, prevI);
curF = initializeMap_withH(prevMap, M, H);
D = single(bwdist(~M));
curF = int32(curF);
numiter = int32(10);
l=0; % special case since this no pyramid is used here.
[outI, outF] = mex_fillOneLevel(curF, I, M, D, l, useLineConstr, numiter); 
figure; imshow(outI);

tel = toc(tst);
fprintf('time of fillOneImage_withinit = %.4f', tel);
% fprintf('level %d end\n', l);

% figure, imshow(resI);
end