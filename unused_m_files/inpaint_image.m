dbstop if error
%% read image and mask
% suppose the image is hxwx3, mask is hxwx1
I = imread('images/halfim.png');
M = imread('images/mask_more_bin.png');
if ~ismatrix(M)
    M = rgb2gray(M);
end
I = im2single(I);
M = im2single(M);
M = M>0;

%% construct pyramid of an image
[pyI, pyM] = constructPyr(I, M); % the levels of the pyramid is computed rather than assigned.

%% inpainting through the pyramid
useLineConstr = 0;
[pyI, F] = fillImagePyr_withmex(pyI, pyM, useLineConstr);

resI = pyI{length(pyI)};

figure, imshow(resI);