function [pyI,F] = fillImagePyr_withmex(pyI, pyM, pynuM, usrLines)
% pyI and pyM are cell vectors containing the image pyramid
% In video inpainting, only the first frame (or the keyframe) needs to be
% processed through this function, others can just be filled using
% fillOneLevel() (see below)
% 一定使用直线，usrLines是用户定义的直线，不用直线的话这个值等于-1。
% 在第一帧之后，就不使用pynuM了，forward的Fhat就直接用其它的能量成分优化。
t1 = tic;
L = length(pyM);
curF = [];
numitertop = floor((linspace(200, 20, L)));
params.alphaSp = 0.025;
params.alphaAp = 0.35;
params.cs_imp = 1;
params.cs_rad = 20;


% linesPyr = cell(L, 1);
linesPyr = pyrLines(pyI, pyM, usrLines, L); % detect lines that are near to or cut the mask

for l = 1:L
    pyI{l} = maskImage(pyI{l}, pyM{l});
    if(l==0)
        curF = initializeMap(curF, pyM{l} );
    else
        curF = vec_initMap(curF, pyM{l});
%         showF(pyI{l}, pyM{l}, curF);
    end
    fprintf('level %d\n', l);
    tic
    D = single(bwdist(~pyM{l}));
    curF = int32(curF);
    numiter = int32(numitertop(l));
%     [pyI{l}, curF] = mex_fillOneLevel( curF, pyI{l}, pyM{l}, D, l, useLineConstr, numiter );
    curF = permute(curF, [3 1 2]); pyI{l} = im2uint8(permute(pyI{l}, [3 1 2]));
    [retI, curF] = mex_fillOneLevel_withline( curF, pyI{l}, pyM{l}, pynuM{l}, D, l, linesPyr{l}, numiter, params );
    curF = ipermute(curF, [3 1 2]); pyI{l} = im2single(ipermute(retI, [3 1 2]));
    %     if(mod(l,2)==1 && l<L)
%         showF(pyI{l}, pyM{l}, curF);
%     end
%     showF(pyI{l}, pyM{l}, curF);
    toc
    fprintf('level %d end\n', l);
end
% showF(pyI{l}, pyM{l}, curF);
toc(t1)
F = curF;
end