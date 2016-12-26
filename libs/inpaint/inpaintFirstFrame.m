function [I, F, C, OF] = inpaintFirstFrame(I)
% inpaint an image _I_ from scratch
OF = object_tracking([], I, 1, [], []); % initialize OF

I = im2single(I);

[M, nouseM, usrLn, C] = getConstraintsGUI(I);

if ~ismatrix(M)
    M = rgb2gray(M);
    nouseM = rgb2gray(nouseM);
end
M = M>0;
nouseM = nouseM>0;

%% construct pyramid of an image
[pyI, pyM, pynuM] = constructPyr(I, M, nouseM); % the levels of the pyramid is computed rather than assigned.

%% inpainting through the pyramid
% useLineConstr;
[pyI, F] = fillImagePyr_withmex(pyI, pyM, pynuM, usrLn);

I = pyI{length(pyI)};

% figure, imshow(resI);

end