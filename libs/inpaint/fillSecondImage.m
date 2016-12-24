function [I, F] = fillSecondImage(I, M, initF, refI)
% Fill image _I_ with hole _M_, initialized mapping _F_ and reference
% model _refI_
% No pyramid processing
% No line constraint

params.alphaSp = 0.025;
params.alphaAp = 0.25;

D = single(bwdist(~M));
F = initF;
numiter = 20;
lines = [];
[I, F] = mex_fillSecondImage(F, I, M, D, refI, lines, numiter, params);

end