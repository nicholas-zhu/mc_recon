% noncart_motion_tv_demo
% preparation
clc;
clear;
addpath(genpath('ESPIRiT'));
addpath(genpath('fminlbfgs_version2c'));
addpath /Users/xuchengzhu/Documents/Graduate/bart/matlab/
load nocart_data.mat;
load smap.mat

% binning, how many bins
Res_d = diff(filter([0.1 0.2 0.4 0.2 0.1],1,Res_Signal));
Res = sign([Res_d(1);Res_d]).*Res_Signal;
motion_signal = [Res,(1:length(Res_Signal))'];
motion_signal = sortrows(motion_signal,1);
motion_signal = motion_signal(:,2);
nbin = 4;

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
kdata_sort = kdata_sort.*repmat(sqrt(w),1,1,size(b1,3),1);
I = reshape(E'*kdata_sort,240,320,1,nbin);
% Smap estimation

%% admm
% argmin_x ||Ex-d||_2 + lambda/2 * ||x||_TVx
% argmin_x ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)
% argmin_z ||Ex-d||_2 + rho/2 * ||x-z||_2 +y'*(x-z)

iter = 0;
% dual update
rho = 1;
lambda = .2; %depends on artifacts level

% E, TVx define
mTVx = TV2d(nbin);
if exist('B','var')
    mTVx = TV2dm(nbin,B,F);
end


x_k0 = zeros(size(I));
z_k0 = zeros(size(I));
y_k0 = zeros(size(I));
d = kdata_sort(:);

I = reshape(E'*kdata_sort,240,320,1,nbin);
sizeI = size(I);
sizeK = size(kdata_sort);
while(iter<=70)
    iter = iter+1;
%     L2 normalization
%     L2fun = @(x)1/2 * x'*(E'*(E*x)) - 2 * d'*E*x + y_k0'*(x - z_k0) + rho/2 * (x - z_k0)'*(x - z_k0);
%     L2fun = @(x)1/2 * x'*(mft2(x,mask_m,smap)) - 2 * d'*mft(x,mask_m,smap)...
%         + y_k0'*(x - z_k0) + rho/2 * (x - z_k0)'*(x - z_k0);
%     [x_k1,fval,exitflag,output,grad] = fminlbfgs(L2fun,x_k0);
%     x_k1 = reshape(x_k1,size(mask_m));

    z_k0 = z_k0(:);
    y_k0 = y_k0(:);
    x_k0 = x_k0(:);
    % L2 opt conjugate gradient descent
    cg_iter = 0;
    cg_iterM = 30;
    x_k = x_k0;
    cg_thd = 0.005;
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
    x_k1 = reshape(x_k1,sizeI);
    y_k0 = reshape(y_k0,sizeI);
    % L1 normalization
    temp_z = mTVx*( x_k1 + y_k0/rho);
    temp_z1 = temp_z(:,:,:,:,2);
    temp_z(:,:,:,:,1) = wthresh(temp_z(:,:,:,:,1),'s',lambda/rho*max(abs(temp_z1(:))));
    z_k1 = mTVx'*temp_z;
    % dual update
    y_k1 = y_k0 + rho * (x_k1 -z_k1);
    
    % all update
    x_k0 = x_k1;
    y_k0 = y_k1;
    z_k0 = z_k1; 
    
end
% state average
writecfl('result_noTV_5',x_k1);
% state sequential
% state loop
%
% registration
%%
I = single(squeeze(abs(x_k0)));
gamma = 1;
I = single(I);
I = I./max(I(:));
I = reshape(imadjust(I(:,:),[0,1],[0,1],gamma),size(I));

Options.Similarity = 'sd';
Options.Spacing = [64 64];
Options.Penalty = 1e-3;
Options.Registration = 'NonRigid';
for i = 1:size(I,3)
    [Ireg(:,:,i),O_trans,Spacing,M,B(:,:,:,i),F(:,:,:,i)] = image_registration(I(:,:,i),I(:,:,1),Options);
end
%%
% TV redesign
% algorithm design


% low rank model plugin
