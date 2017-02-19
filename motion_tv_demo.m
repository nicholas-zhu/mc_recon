% preparation
clc;
clear;
addpath(genpath('ESPIRiT'));
addpath(genpath('fminlbfgs_version2c'));

load recon_data.mat
% d = repmat(DATAAc,1,1,1,1,4).*repmat(mask_m,1,1,1,6,1);
% mask_m = (mask_m(:,:,:,:,[1,4])+mask_m(:,:,:,:,[2,3]));
d = repmat(DATAAc,1,1,1,1,4).*repmat(mask_m,1,1,1,6,1);
% initialization x, y, 
% data normalization

%% admm
% argmin_x ||Ex-d||_2 + lambda/2 * ||x||_TVx
% argmin_x ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)
% argmin_z ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)

iter = 0;
% dual update
rho = 1;
lambda = 0.002;

% E, TVx define
E = FFTx(smap,mask_m);
mTVx = TV();

x_k0 = zeros(size(mask_m));
z_k0 = zeros(size(mask_m));
y_k0 = zeros(size(mask_m));
d = d(:);

while(iter<=20)
    iter = iter+1;
%     L2 normalization
    z_k0 = z_k0(:);
    y_k0 = y_k0(:);
    x_k0 = x_k0(:);
    % L2 opt conjugate gradient descent
    cg_iter = 0;
    cg_iterM = 100;
    x_k = x_k0;
    cg_thd = 0.01;
    beta_k = 0;
    grad_k = (E'*(E*x_k0-d)) + rho*(x_k0-z_k0) + y_k0;
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
        grad_k1 = (E'*(E*x_k1-d)) + rho*(x_k1-z_k0) + y_k0;
        %beta_k = d_k'*(E'*(E*grad_k1))/(d_k'*(E'*(E*d_k))+eps);
        beta_k = grad_k1'*grad_k1/(grad_k'*grad_k+eps);
        
        d_k = -grad_k1+beta_k*d_k;
        grad_k = grad_k1;
        x_k = x_k1;
    end
    all_value = norm(x_k1(:)-x_k0(:))/(norm(x_k1(:))+eps);
    if (all_value<0.005)
        break;
    end
    fprintf(['cost function value:',num2str(cg_value),'\n']);
    x_k1 = reshape(x_k1,size(mask_m));
    y_k0 = reshape(y_k0,size(mask_m));
    % L1 normalization
    temp_z = mTVx*( x_k1 + y_k0/rho );
    temp_z1 = wthresh(temp_z(:,:,:,:,:,1),'s',lambda/rho*max(abs(temp_z(:))));
    z_k1 = temp_z1+temp_z(:,:,:,:,:,2);
    
    % dual update
    y_k1 = y_k0 + rho * (x_k1 -z_k1);
    
    % all update
    x_k0 = x_k1;
    y_k0 = y_k1;
    z_k0 = z_k1;    
end
save('resultc.mat','x_k1','z_k1');

