function res =LR2dn(bsize)
% improved TV constraint

res.adjoint = 0;
res.bsize = bsize;
res = class(res,'LR2dn');