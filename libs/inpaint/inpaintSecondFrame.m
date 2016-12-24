function [I, F, C, H, OF] = inpaintSecondFrame(I, refI, refF, C, OF)
% Inpaint image _I_ with a refence fram _refI_, between which there exist a
% homography _H_ which satisfies p = H*p1, p \in _I_ and p1 in _refI_. 
% _refF_ is the of _refI_, _C_ the contour points of the object to be
% removed in _refI_, _OF_ the optical flow which is used in tracking.

%% Pipeline: 
% 1. tracking and determine H
% 2. mapping forwarding to get hatF
% 3. refine hatF for I to get F and I (with mex file)

%% Pt. 1

%% Pt. 2

%% Pt. 3

end