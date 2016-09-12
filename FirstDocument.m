% The document..

%% The most basic steps
% First, write an image inpainting algorithm based on three-fold energy
% function. (without parallel)

% The image inpainting algorithm is based on pyramid. To propagate the
% mapping from a coarser level to a finer level, I think I may just use an
% easy interpolation method.

% One HARD step is the optimization of the energy function which involves a
% random sampling step and a propagation step. The original authors used a
% parallel framework to do this but I think I will just do it sequentially
% for the whole inpainting area.

%% See what it will bring us without constraints and do the following
% The paper used only straight line constraint so I will follow them. That
% should be based on Hough transformation.

% Then, add the homography estimation and the reference model to
% incorporate video inpainting. 
% Since we do not add