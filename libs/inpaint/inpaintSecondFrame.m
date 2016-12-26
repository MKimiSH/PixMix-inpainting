function [I, M, F, C, H, OF, OFOBJ] = inpaintSecondFrame(I, lastI, lastF, lastM, lastC, OF, OFOBJ)
% Inpaint image _I_ with a refence fram _lastI_, between which there exist a
% homography _H_ which satisfies p = p1*H, p \in _I_ and p1 in _lastI_ (row vectors). 
% _lastF_ is the px mapping of _lastI_, _C_ the contour points of the object to be
% removed in _lastI_, _OF_ the optical flow which is used in tracking.

%% Pipeline: 
% 1. Track and determine H (BE careful of H'!!)
% 2. Mapping forwarding to get initF, lightness adjustment
% 3. Refine initF to get F and I, and (if necessary) blend the edges

%% Pt. 1
tic;
[H, C, M, OFOBJ, ~, ~, OF] = object_tracking_novec(lastI, I, 0, lastC, OFOBJ);
t1 = toc; fprintf('%.4f sec for tracking\n', t1);
%% Pt. 2
H = H';
tic;
% [M] = getMFromH(lastM, H);
[refI, initF] = forwardF(I, M, lastI, lastF, lastM, H);
% figure, imshow(M);
t2 = toc; fprintf('%.4f sec for forwarding\n', t2);
% refI = lightnessAdjust(refI, M, C, lastI, lastC);
%% Pt. 3
tic;
[I, F] = fillSecondImage(I, M, initF, refI);
close;
figure, imshow(I);
t3 = toc; fprintf('%.4f sec for inpainting\n', t3);
end

function [M] = getMFromH(lastM, H)

[row, col] = find(lastM);
row = row(1:2:end); col = col(1:2:end);
P = [col, row, ones(length(row), 1)];

P = P*H;

rrow = round(P(:,2) ./ P(:,3));
rcol = round(P(:,1) ./ P(:,3));

M = zeros(size(lastM), 'logical');

M(sub2ind(size(M), rrow, rcol)) = 1;
M = imclose(M, strel('disk', 3));

end