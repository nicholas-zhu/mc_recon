function res =TV2dm(L,B,F)
% improved TV constraint

res.F = F;
res.B = B;
res.adjoint = 0;
res.w = ones(L,1);
res = class(res,'TV2dm');
