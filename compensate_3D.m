%% motion compensated recon
addpath(genpath('../code'));

%% preparation
% data structure
ksp = readcfl('../Navigator/ksp');
smap = readcfl('../Navigator/smap');
smap = smap(:,:,:,:,1);
Img = readcfl('../Navigator/pics_img');
Img = Img(:,:,:,:,1,:);
ksp = squeeze(permute(ksp,[2 3 1 4 5 6]));
smap = squeeze(permute(smap,[2 3 1 4]));
Img = permute(Img,[2 3 1 4 6 5]);

% x,y,z(full sample),coil, motion
%% motion field estimation
% downsampling 4= readc
B = readcfl('B');
F = readcfl('F');
isreg = 0;
nbin = size(Img,5);
if(isreg)
downsampling = 2;
Img = abs(Img)./max(abs(Img(:)))+eps;
nbin = size(Img,5);
sizeH = size(Img(:,:,:,1));
sizeL = floor(sizeH/downsampling);
[vy,vx,vz] = ndgrid(linspace(1,sizeH(1),sizeL(1)),...
    linspace(1,sizeH(2),sizeL(2)),...
    linspace(1,sizeH(3),sizeL(3)));
Img_L = zeros([sizeL,1,nbin]);
for i = 1:nbin
    Img_L(:,:,:,1,i) = interp3(Img(:,:,:,i),vx,vy,vz,'cubic');
end
Img_L = Img_L.^0.5;
%%
ref_N = 4;
Options.Similarity = 'sd';
Options.Spacing = [16,8,16];
Options.Penalty = 10e-5;
Options.Registration = 'NonRigid';
Ireg = zeros(size(Img_L));

B_L = zeros([size(Img_L(:,:,:,1)),3,nbin]);
F_L = zeros([size(Img_L(:,:,:,1)),3,nbin]);
% regi
for i = 1:nbin
    [Ireg(:,:,:,1,i),~,~,~,B_L(:,:,:,:,i),F_L(:,:,:,:,i)] = image_registration(Img_L(:,:,:,1,i),Img_L(:,:,:,1,ref_N),Options);
end
writecfl('Ireg',Ireg);
%%
B = zeros([sizeH,3,nbin]);
F = zeros([sizeH,3,nbin]);
[vy,vx,vz] = ndgrid(linspace(1,sizeL(1),sizeH(1)),...
    linspace(1,sizeL(2),sizeH(2)),...
    linspace(1,sizeL(3),sizeH(3)));

for i = 1:nbin
    for j = 1:3
        B(:,:,:,j,i) = downsampling*interp3(B_L(:,:,:,j,i),vx,vy,vz,'cubic');
        F(:,:,:,j,i) = downsampling*interp3(F_L(:,:,:,j,i),vx,vy,vz,'cubic');
    end
end
writecfl('B',B);
writecfl('F',F);
end
% interpolate 4
%% reconstruction framework
% motion compensated 
mask = abs(ksp(:,:,:,1,:))>0;
E = FFTx(smap,mask);
mTVx = TVm(nbin,B,F);
Wav = Wavelet('Daubechies_TI',8,8);
mTVs = TVs([1 2]);
sizeI2 = size(Img);
I = reshape(E'*ksp,sizeI2);
d = ksp(:);

rho = 1;
rho1 = 1;
A = @(x,flag)(E'*(E*x)+(rho+rho1)*x);
Ad = E'*d;
lambda1 = 0.05*max(abs(Ad(:)));
lambda2 = 0.01;

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
iter = 0;
tic;
while(iter<=30)
    iter = iter+1;
    % L2 opt conjugate gradient descent
    cg_iterM = 30;
    tol = 1e-3;
    b = Ad + rho *(z_k0(:)-u_k0(:)) + rho1 *(z1_k0(:)-u1_k0(:));
    x_k1 = lsqr(A,b,tol,cg_iterM,[],[],x_k0(:));
    x_k1 = reshape(x_k1,sizeI2);
    
    % L1 normalization
    temp_z = mTVx*( x_k1+ u_k0);
    temp_z1 = temp_z(:,:,:,:,:,1);
    temp_z(:,:,:,:,:,1) = wthresh(temp_z1,'s',lambda1);
    z_k1 = mTVx'*temp_z;
    
    % L1 normalization
%     temp_z = TVs*(x_k1+ u1_k0);
%     temp_z1 = temp_z(:,:,:,:,:,1,:);
%     temp_z(:,:,:,:,:,1,:) = wthresh(temp_z1,'s',lambda1);
    temp1 = x_k1 + u1_k0;
    temp1 = temp1(:,:,:);
    tempZ = [256,256,size(temp1,3)];
    tempS = size(temp1);
    temp1 = zpad(temp1,tempZ);
    temp1 = Wav*(temp1);
    temp1 = Wav'*(wthresh(temp1,'s',lambda2));
    temp1 = crop(temp1,tempS);
    z1_k1 = reshape(temp1,sizeI2);
%    temp_z = mTVs'*temp_z;
%    z1_k1 = reshape(temp_z,sizeI2);
    % dual update
    u_k1 = u_k0 + (x_k1 -z_k1);
    u1_k1 = u1_k0 + (x_k1 -z1_k1);
    
    % all update
    x_k0 = x_k1;
    u_k0 = u_k1;
    z_k0 = z_k1; 
    u1_k0 = u1_k1;
    z1_k0 = z1_k1;
    writecfl('temp_3D_mc1',x_k0);
end
% L1 wavelet
% data consistancy
toc
