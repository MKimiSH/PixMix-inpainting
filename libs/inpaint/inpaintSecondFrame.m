function [I, F, C, H, OF] = inpaintSecondFrame(I, lastI, lastF, lastM, lastC, OF)
% Inpaint image _I_ with a refence fram _lastI_, between which there exist a
% homography _H_ which satisfies p = p1*H, p \in _I_ and p1 in _lastI_ (row vectors). 
% _lastF_ is the px mapping of _lastI_, _C_ the contour points of the object to be
% removed in _lastI_, _OF_ the optical flow which is used in tracking.

%% Pipeline: 
% 1. Track and determine H (BE careful of H'!!)
% 2. Mapping forwarding to get initF, lightness adjustment
% 3. Refine initF to get F and I, and (if necessary) blend the edges

%% Pt. 1
[H, C, M, OF] = object_tracking(lastI, I, 0, lastC, OF);
%% Pt. 2
[refI, initF] = forwardF(I, M, lastI, lastF, lastM, H);
refI = lightnessAdjust(refI, C, lastI, lastC);
%% Pt. 3
[I, F] = fillSecondImage(I, M, initF, refI);

end