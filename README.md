## PixMix-inpainting
A MATLAB reimplementation of TVCG14 paper High-Quality Real-Time Video Inpainting with PixMix (J. Herling and W. Broll)

### Prerequisites
- MATLAB with Optical Flow and Tracking
- mex compiler

### Usage
0. add path recursively to MATLAB with `pathtool`
1. compile mex files

```
mex libs\inpaint\mex_fillOneLevel_withline.cpp
mex libs\inpaint\mex_fillSecondImage.cpp
mex libs\detect\mex_detect_smallst.cpp
```

2. run demo.m->没写说个蛋。


### Before 12.27
I won't incorporate the parallel machanism in the paper, so the speed won't reach real-time. 

Now inpainting of 1st frame ~3s (1280x720 with ~30000 pixels missing) and next frames ~1s (little difference from 1st frame). Bottleneck is tracking, which has a lot to improve

