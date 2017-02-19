function res =TV2d(L)
% improved TV constraint

res.adjoint = 0;
res.w = hamming(L);
res = class(res,'TV2d');