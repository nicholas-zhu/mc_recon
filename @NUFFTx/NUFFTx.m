function res = NUFFTx(k,w,S)

% Multi-bin 2D NUFFT
% based on NUFFT toolbox from Jeffery Fessler
% Inputs:
%  k: fe,pe,1,motion
%  w: fe,pe,1,motion
%  

assert(sum((size(k)~=size(w)))==0,'Traj, Weight not match!');
Nd = size(S(:,:,1));
Jd = [6,6];
Kd = floor(Nd*1.5);
n_shift = Nd/2;

for i = 1:size(k,4)
    km = k(:,:,1,i);
    om = [real(km(:)),imag(km(:))]*2*pi;
    res.st{i} = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser');
end

res.adjoint = 0;
res.imSize = size(S(:,:,1));
res.imSize2 = [res.imSize,1,size(k,4)];
res.dataSize = size(k);
res.dataSize(3) = size(S,3);
res.w = w;
res.S = S;
res = class(res,'NUFFTx');


    
    