function ratio = non_cart_artifacts_level(traj,b1)

w_1 = ones(size(traj));
sizeI = [size(b1,1),size(b1,2)];
Ep = NUFFTx(traj,w_1,double(b1));
I = zeros(size(b1,1),size(b1,2));
I(floor(sizeI/2)+1) = 1;
ksp = Ep*I(:);

cg_iterM = 30;
tol = 2e-2;
b = Ep'*ksp;
Ex = @(x,flag)(Ep'*(Ep*x));
I_1 = lsqr(Ex,b,tol,cg_iterM,[],[],I(:)*0);
ratio = max((I(:)==0).*abs(I_1(:)-I(:)));


