function pyI = fillImagePyr(pyI, pyM, useLineConstr)
% pyI and pyM are cell vectors containing the image pyramid
% In video inpainting, only the first frame (or the keyframe) needs to be
% processed through this function, others can just be filled using
% fillOneLevel() (see below)

L = length(pyM);
curF = [];

for l = 1:L
    curF = initializeMap(curF, pyM{l} );
    fprintf('level %d', l);
    tic
    [pyI{l}, curF] = fillOneLevel( curF, pyI{l}, pyM{l}, l, useLineConstr ); 
    imshow(pyI{l});
    toc
%     prevF = curF;
end

end