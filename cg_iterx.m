% noncart_motion_tv_demo
% preparation
clc;
clear;
addpath(genpath('ESPIRiT'));
addpath(genpath('fminlbfgs_version2c'));
load nocart_data.mat;

% binning, how many bins
motion_signal = [Res_Signal,(1:length(Res_Signal))'];
motion_signal = sortrows(motion_signal,1);
motion_signal = motion_signal(:,2);
nbin = 10;

k_sort = double(k(:,motion_signal));
kdata_sort = double(kdata(:,motion_signal,:));


k_sort = reshape(k_sort,size(k,1),size(k,2)/nbin,1,nbin);
kdata_sort = permute(reshape(kdata_sort,size(k,1),...
    size(k,2)/nbin,nbin,size(kdata,3)),[1 2 4 3]);
% N1 = 200/256;
% N2 = 400/256;
N1 = 1;
N2 = 1;
% w = sqrt(abs(k_sort));
w = ones(size(k_sort));
k_sort = real(k_sort)*N1+1i*imag(k_sort)*N2;
E = NUFFTx(k_sort,w,double(b1));
kdata_sort = kdata_sort.*repmat(sqrt(w),1,1,size(b1,3),1);
I = reshape(E'*kdata_sort,256,256,1,nbin);
% Smap estimation

%% admm
% argmin_x ||Ex-d||_2 + lambda/2 * ||x||_TVx
% argmin_x ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)
% argmin_z ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)

iter = 0;
% dual update
rho = 0.5;
lambda = 0.02;

% E, TVx define
mTVx = TV2d();

x_k0 = zeros(size(I));
z_k0 = zeros(size(I));
y_k0 = zeros(size(I));
d = kdata_sort(:);

sizeI = size(I);
sizeK = size(kdata_sort);
while(iter<=1)
    iter = iter+1;
    z_k0 = z_k0(:);
    y_k0 = y_k0(:);
    x_k0 = x_k0(:);
    % L2 opt conjugate gradient descent
    cg_iter = 0;
    cg_iterM = 100;
    x_k = x_k0;
    cg_thd = 0.01;
    beta_k = 0;
    grad_k = (E'*(E*x_k0-d)) + rho*(x_k0);
    d_k = -grad_k;
    while(cg_iter<cg_iterM)
        cg_iter = cg_iter+1;
        lambda_k = -d_k'*grad_k/(d_k'*(E'*(E*d_k)+rho*(d_k'*d_k))+eps);
        x_k1 = x_k + lambda_k * d_k;
        
        cg_value = norm(x_k1-x_k)/(norm(x_k1)+eps);
        fprintf(['cost function value:',num2str(cg_value),'\n']);
        if cg_value<cg_thd
            break;
        end
        grad_k1 = (E'*(E*x_k1-d)) + rho*(x_k1);
        %beta_k = d_k'*(E'*(E*grad_k1))/(d_k'*(E'*(E*d_k))+eps);
        beta_k = (grad_k1'*grad_k1)/(grad_k'*grad_k+eps);
        
        d_k = -grad_k1+beta_k*d_k;
        grad_k = grad_k1;
        x_k = x_k1;
    end
    all_value = norm(x_k1(:)-x_k0(:))/(norm(x_k1(:))+eps);
    if (all_value<0.005)
        break;
    end
    fprintf(['cost function value:',num2str(cg_value),'\n']);
    x_k1 = reshape(x_k1,sizeI);
    
    % all update
    x_k0 = x_k1;   
end
% NUFFT design
% TV redesign
% algorithm design
% low rank model plugin
