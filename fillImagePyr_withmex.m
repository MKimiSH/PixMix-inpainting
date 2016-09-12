function [pyI,F] = fillImagePyr_withmex(pyI, pyM, useLineConstr)
% pyI and pyM are cell vectors containing the image pyramid
% In video inpainting, only the first frame (or the keyframe) needs to be
% processed through this function, others can just be filled using
% fillOneLevel() (see below)

L = length(pyM);
curF = [];

for l = 1:L
    curF = initializeMap(curF, pyM{l} );
    fprintf('level %d\n', l);
    tic
%     [R,C] = size(pyM{l});
    D = single(bwdist(~pyM{l}));
    curF = int32(curF);
    numiter = int32(10);
    [pyI{l}, curF] = mex_fillOneLevel( curF, pyI{l}, pyM{l}, D, l, useLineConstr, numiter ); 
%     imshow(pyI{l});
    toc
    fprintf('level %d end\n', l);
%     prevF = curF;
end
F = curF;
end