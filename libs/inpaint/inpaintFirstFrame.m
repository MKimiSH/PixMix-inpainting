function [I, M, F, C, OF, OFobj] = inpaintFirstFrame(I, M, nouseM, usrLn, C)
% inpaint an image _I_ from scratch
[~, ~, ~, OFobj, ~, ~, OF] = object_tracking_novec([], I, 1, C, []); % initialize OF

I = im2single(I);

% [M, nouseM, usrLn, C] = getConstraintsGUI(I);

if ~ismatrix(M)
    M = rgb2gray(M);
    nouseM = rgb2gray(nouseM);
end
M = M>0;
nouseM = nouseM>0;

M = imdilate(M,strel('disk',5));
%% construct pyramid of an image
[pyI, pyM, pynuM] = constructPyr(I, M, nouseM); % the levels of the pyramid is computed rather than assigned.

%% inpainting through the pyramid
% useLineConstr;
[pyI, F] = fillImagePyr_withmex(pyI, pyM, pynuM, usrLn);

I = pyI{length(pyI)};

figure, imshow(I);
imwrite(I, 'out1.jpg');

end