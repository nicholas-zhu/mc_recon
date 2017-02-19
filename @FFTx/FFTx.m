function res = FFTx(Smap, Mask)

% Multi-bin 3D FFT
% Inputs:
%  
%  Smap: kx, ky, kz, S, 1
%  Mask: kx, ky, kz, 1, M
res.adjoint = 0;
res.Smap = Smap;
res.sizeS = size(Smap,4);
res.nMask = size(Mask,5);
res.Mask = Mask;
res.sizeM = size(Mask,5);
res = class(res,'FFTx');


