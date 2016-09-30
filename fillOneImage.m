function resI = fillOneImage(I, M, useLineConstr)
if ~ismatrix(M)
    M = rgb2gray(M);
end
I = im2single(I);
M = im2single(M);
M = M>0;

%% construct pyramid of an image
[pyI, pyM] = constructPyr(I, M); % the levels of the pyramid is computed rather than assigned.

%% inpainting through the pyramid
% useLineConstr;
pyI = fillImagePyr_withmex(pyI, pyM, useLineConstr);

resI = pyI{length(pyI)};

% figure, imshow(resI);
end