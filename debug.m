% noncart_motion_tv_demo
% preparation
tic;
clc;
clear;
addpath(genpath('/Users/xuchengzhu/Documents/Graduate/Research/Projects/Motion_Recon/reg_recon'));
addpath /Users/xuchengzhu/Documents/Graduate/bart/matlab/
load nocart_dce_data.mat;
load smap_dce.mat;

motion_signal = [(1:length(Res_Signal))',Res_Signal];
motion_signal = sortrows(motion_signal,2);
motion_signal = motion_signal(:,1);


% Res_Signal = 1:length(Res_Signal);
% motion_signal = Res_Signal;
% motion_signal = sortrows(motion_signal,1);

nbin = 5;

k_sort = double(k(:,motion_signal));
kdata_sort = double(kdata(:,motion_signal,:));
b1 = smap;


k_sort = reshape(k_sort,size(k,1),size(k,2)/nbin,1,nbin);
kdata_sort = permute(reshape(kdata_sort,size(k,1),...
    size(k,2)/nbin,nbin,size(kdata,3)),[1 2 4 3]);
N1 = 1;
N2 = 1;

% w = sqrt(abs(k_sort));
w = ones(size(k_sort));
k_sort = real(k_sort)/N1+1i*imag(k_sort)/N2;
E = NUFFTx(k_sort,w,double(b1));
% kdata_sort = kdata_sort.*repmat(sqrt(w),1,1,size(b1,3),1);
I = reshape(E'*kdata_sort,240,320,1,nbin);
% Smap estimation

%% admm
% argmin_x ||Ex-d||_2 + lambda/2 * ||x||_TVx
% argmin_x ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)
% argmin_z ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)

iter = 0;
% dual update
rho = 1;
lambda = .1; %depends on artifacts level
w = repmat(abs(k_sort),1,1,size(b1,3),1);
w = w(:);
d = kdata_sort(:);

% E, TVx define
mTVx = TV2do;
A = @(x,flag)(E'*(w.*(E*x))+(rho)*x);
Ad = E'*(w.*d);

if exist('B','var')
    mTVx = TV2dm(nbin,B,F);
end

if ~exist('x_k0','var')
    x_k0 = zeros(size(I));
    z_k0 = zeros(size(I));
    u_k0 = zeros(size(I));
end
I = reshape(E'*kdata_sort,240,320,1,nbin);
sizeI1 = size(I);
sizeK = size(kdata_sort);
while(iter<=20)
    iter = iter+1;
    z_k0 = z_k0(:);
    u_k0 = u_k0(:);
    x_k0 = x_k0(:);
    
    % L2 opt conjugate gradient descent
    cg_iterM = 20;
    x_k = x_k0;
    tol = 1e-4;
    b = Ad + rho *(z_k0-u_k0);
    x_k1 = lsqr(A,b,tol,cg_iterM,[],[],x_k0);
%     while(cg_iter<cg_iterM)
%         cg_iter = cg_iter+1;
%         lambda_k = -d_k'*grad_k/(d_k'*(E'*(w.*(E*d_k))+rho*(d_k'*d_k))+eps);
%         x_k1 = x_k + lambda_k * d_k;
%         
%         cg_value = norm(x_k1-x_k)/(norm(x_k1)+eps);
%         fprintf(['cost function value:',num2str(cg_value),'\n']);
%         if cg_value<cg_thd
%             break;
%         end
%         grad_k1 = (E'*(w.*(E*x_k1-d))) + rho*(x_k1-z_k0+u_k0);
%         %beta_k = d_k'*(E'*(E*grad_k1))/(d_k'*(E'*(E*d_k))+eps);
%         beta_k = grad_k1'*grad_k1/(grad_k'*grad_k+eps);
%         
%         d_k = -grad_k1+beta_k*d_k;
%         grad_k = grad_k1;
%         x_k = x_k1;
%     end
    x_k1 = reshape(x_k1,sizeI1);
    u_k0 = reshape(u_k0,sizeI1);
    % L1 normalization
    temp_z = mTVx*( x_k1 + u_k0);
    temp_z1 = temp_z(:,:,:,:,2);
    temp_z(:,:,:,:,1) = wthresh(temp_z(:,:,:,:,1),'s',lambda/rho*max(abs(temp_z1(:))));
    z_k1 = mTVx'*temp_z;
    
    % dual update
    u_k1 = u_k0 + (x_k1 -z_k1);
    
    % all update
    x_k0 = x_k1;
    u_k0 = u_k1;
    z_k0 = z_k1; 
    
end
writecfl('mc_dce',x_k1);