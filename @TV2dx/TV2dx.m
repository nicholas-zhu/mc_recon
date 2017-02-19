function res =TV2dx(B,F)
% improved TV constraint
% no weighting
res.adjoint = 0;
res.B = B;
res.F = F;
res = class(res,'TV2dx');