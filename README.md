## PixMix-inpainting
A MATLAB reimplementation of TVCG14 paper High-Quality Real-Time Video Inpainting with PixMix

I won't incorporate the parallel machanism in the paper, so the speed won't reach real-time, but it's still quite fast, about several seconds per frame. 

mex libs\inpaint\mex_fillSecondImage.cpp
mex libs\inpaint\mex_fillOneLevel_withline.cpp