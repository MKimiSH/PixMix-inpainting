function [I, F] = fillSecondImage(I, M, initF, refI)
% Fill image _I_ with hole _M_, initialized mapping _F_ and reference
% model _refI_
% No pyramid processing
% No line constraint

params.alphaSp = 0.025;
params.alphaAp = 0.95;

D = single(bwdist(~M));
numiter = int32(20);
lines = [];
F = permute(int32(initF), [3 1 2]);
I = permute(im2uint8(I), [3 1 2]);
prefI = permute(im2uint8(refI), [3 1 2]);
[I, F] = mex_fillSecondImage(F, I, M, D, prefI, lines, numiter, params);
F = ipermute(F, [3 1 2]);
I = im2single(ipermute(I, [3 1 2]));

end