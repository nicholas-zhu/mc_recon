function res =TVm(L,B,F)
%  3DTV motion constraint

res.F = F;
res.B = B;
res.adjoint = 0;
res.w = hamming(L);
res = class(res,'TVm');