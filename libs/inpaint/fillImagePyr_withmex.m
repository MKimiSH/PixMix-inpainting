function [pyI,F] = fillImagePyr_withmex(pyI, pyM, useLineConstr)
% pyI and pyM are cell vectors containing the image pyramid
% In video inpainting, only the first frame (or the keyframe) needs to be
% processed through this function, others can just be filled using
% fillOneLevel() (see below)
t1 = tic;
L = length(pyM);
curF = [];
numitertop = floor((linspace(200, 10, L)));
params.alphaSp = 0.025;
params.alphaAp = 0.25;
params.cs_imp = 1;
params.cs_rad = 20;

linesPyr = cell(L, 1);
if useLineConstr>0
    linesPyr = pyrLines(pyI, pyM, L); % detect lines that are near to or cut the mask
end

for l = 1:L
    pyI{l} = maskImage(pyI{l}, pyM{l});
    curF = initializeMap(curF, pyM{l} );
    if(l==1)
%         showF(pyI{l}, pyM{l}, curF);
    end
    fprintf('level %d\n', l);
    tic
%     [R,C] = size(pyM{l});
    D = single(bwdist(~pyM{l}));
    curF = int32(curF);
    numiter = int32(numitertop(l));
%     [pyI{l}, curF] = mex_fillOneLevel( curF, pyI{l}, pyM{l}, D, l, useLineConstr, numiter );
    [pyI{l}, curF] = mex_fillOneLevel_withline( curF, pyI{l}, pyM{l}, D, l, linesPyr{l}, numiter, params );
    if(mod(l,2)==1 && l<L)
        showF(pyI{l}, pyM{l}, curF);
    end
    toc
    fprintf('level %d end\n', l);
%     prevF = curF;
end
showF(pyI{l}, pyM{l}, curF);
toc(t1)
F = curF;
end