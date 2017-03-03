% noncart_motion_tv_demo
% preparation
tic;
clc;
clear;
addpath(genpath('/working/larson/xzhu/code'));
load nocart_dce_data.mat;
load smap_dce.mat;
smap = squeeze(readcfl('smap_dce'));
motion_signal = [(1:length(Res_Signal))',Res_Signal];
motion_signal = sortrows(motion_signal,2);
motion_signal = motion_signal(:,1);

nbin = 5;
mc_flag = 1;

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
w = repmat(abs(k_sort),1,1,size(b1,3),1);
w = w(:);
d = kdata_sort(:);

% E, TVx define
Ad = E'*(w.*d);

I = reshape(E'*kdata_sort,240,320,1,nbin);
sizeI1 = size(I);
sizeK = size(kdata_sort);
x_k1 = reshape(Ad,sizeI1);

writecfl('mc_dce',x_k1);
%%
if (mc_flag)
I = single(squeeze(abs(x_k1)));
gamma = 0.5;
I = single(I);
I = I./max(I(:));
I = reshape(imadjust(I(:,:),[0,1],[0,1],gamma),size(I));

Options.Similarity = 'sd';
Options.Spacing = [32 32];
Options.Penalty = 1e-3;
Options.MaskMoving = (smap(:,:,1)>eps);
for i = 1:size(I,3)
    [Ireg(:,:,i),O_trans,Spacing,M,B(:,:,:,i),F(:,:,:,i)] = image_registration(I(:,:,i),I(:,:,1),Options);
end

writecfl('temp',Ireg);
end
%% DCE recon
% low rank
% regroup
motion_signal = [];
w_t = [];
dce_state = 11;
Nbin = nbin*dce_state;
for i = 0:(dce_state-1)
    index = i*100+(1:100);
    motion_t = [index',Res_Signal(index)];
    motion_t = sortrows(motion_t,2);
    motion_signal = [motion_signal;motion_t(:,1)];
    w_t = [w_t;(motion_t(:,1)-100*i)];
end
w_t = w_t-51;
w_t1 = linspace(100,20,256);
w_t2 = [w_t1,flip(w_t1,2)]';
w_time = zeros(512,length(w_t));
for i =1:length(w_t)
    w_time(:,i) = exp(-w_t(i).^2/(2*w_t2.^2));
end
w_time = repmat(reshape(w_time,size(k,1),size(k,2)/Nbin,1,Nbin),...
    1,1,size(kdata,3),1);
w_time = w_time(:);
k_sort = double(k(:,motion_signal));
kdata_sort = double(kdata(:,motion_signal,:));
b1 = smap;

Nbin = nbin*dce_state;
sizeI1 = [240,320,1,Nbin];
sizeI2 = [240,320,1,nbin,dce_state];
sizeK1 = [size(k,1),size(k,2)/Nbin,1,Nbin];
sizeK2 = [size(k,1),size(k,2)/Nbin,1,nbin,dce_state];
k_sort = reshape(k_sort,sizeK1);
kdata_sort = permute(reshape(kdata_sort,size(k,1),...
    size(k,2)/Nbin,Nbin,size(kdata,3)),[1 2 4 3]);

iter = 0;
rho = 1;
rho1 = 1;


% E, TVx define
%Bt = repmat(B,1,1,1,dce_state);
%Ft = repmat(F,1,1,1,dce_state);
%mTVx = TV2dm(Nbin,Bt,Ft);
%mLRx = LR2d([6,6],Bt,Ft);
%mTVx defined before
if(mc_flag)
mTVx = TV2dm(nbin,B,F);
else
mTVx = TV2do();
end
mTVx2 = TV2do();
mLRx = LR2dn([4,4]);
w = ones(size(k_sort));
E = NUFFTx(k_sort,w,double(b1));
I = reshape(E'*kdata_sort,240,320,1,Nbin);
w = repmat(abs(k_sort),1,1,size(b1,3),1);
w = 2*w(:);
d = kdata_sort(:);
A = @(x,flag)(E'*(w.*(E*x))+(rho+rho1)*x);
Ad = E'*(w.*d);
writecfl('mc_dce1',reshape(Ad,sizeI2));
smooth = 0.02;
w = (w+smooth)/(1+smooth);
lambda = 0.1*sum(abs(Ad(:)))./sum(abs(Ad(:))>0);
lambda2 = lambda/4;

% init
x_k0 = reshape(Ad,sizeI2);
z_k0 = zeros(sizeI2);
u_k0 = zeros(sizeI2);
z1_k0 = zeros(sizeI2);
u1_k0 = zeros(sizeI2);

x_k1 = zeros(sizeI2);
z_k1 = zeros(sizeI2);
u_k1 = zeros(sizeI2);
z1_k1 = zeros(sizeI2);
u1_k1 = zeros(sizeI2);
%%
while(iter<=30)
    iter = iter+1;
    % L2 opt conjugate gradient descent
    cg_iterM = 30;
    tol = 2e-2;
    b = Ad + rho *(z_k0(:)-u_k0(:)) + rho1 *(z1_k0(:)-u1_k0(:));
    x_k1 = lsqr(A,b,tol,cg_iterM,[],[],x_k0(:));
    x_k1 = reshape(x_k1,sizeI2);
    
    % L1 normalization
    for i = 1:dce_state
        temp_z = mTVx*( x_k1(:,:,:,:,i) + u_k0(:,:,:,:,i));
        temp_z1 = temp_z(:,:,:,:,1);
        temp_z(:,:,:,:,1) = wthresh(temp_z(:,:,:,:,1),'s',lambda);
        z_k1(:,:,:,:,i) = mTVx'*temp_z;
    end
    
    % L* normalization
    for j = 1:nbin
        temp_t = mTVx2*(permute((x_k1(:,:,:,j,:) + u1_k0(:,:,:,j,:)),[1 2 3 5 4]));
        temp_t1 = temp_t(:,:,:,:,1);
        temp_t(:,:,:,:,1) = wthresh(temp_t(:,:,:,:,1),'s',lambda2);
        z1_k1(:,:,:,j,:) = mTVx2'*temp_t;
    end
    
    % dual update
    u_k1 = u_k0 + (x_k1 -z_k1);
    u1_k1 = u1_k0 + (x_k1 -z1_k1);
    
    % all update
    x_k0 = x_k1;
    u_k0 = u_k1;
    z_k0 = z_k1; 
    u1_k0 = u1_k1;
    z1_k0 = z1_k1;
    writecfl('mc_dcekk3t',x_k0); 
end

writecfl('mc_dcekk3',x_k0);
time=toc;
time/60
