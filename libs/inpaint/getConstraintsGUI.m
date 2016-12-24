function [M, nuM, userLines, C] = getConstraintsGUI(I)
% Open a window on which the user can select object that is going to be deleted
% (with a contour _C_ and a mask _M_) and also lines that may help
% inpainting (_userLines_), as well as pixel positions that should not be
% used in the inpainting process for image _I_ (_nuM_)

M = [];
nuM = [];
userLines = [];
C = [];

end