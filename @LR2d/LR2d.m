function res =LR2d(bsize,B,F)
% improved TV constraint

res.adjoint = 0;
res.bsize = bsize;
res.B = B;
res.F = F;
res = class(res,'LR2d');