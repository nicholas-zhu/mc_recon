function res =TV2d(L)
% improved TV constraint

res.adjoint = 0;
%res.w = hamming(L);
res.w = ones(L,1);
res = class(res,'TV2d');
