clc;
clear;
addpath(genpath('ESPIRiT'));
addpath(genpath('fminlbfgs_version2c'));
load nocart_dce_data.mat;
scale = 1;
k = double(k);
kdata = double(kdata);
k = k(247:266,1:150);
kdata = kdata(247:266,1:150,:);

w = abs(k);
k = k/scale;
smap_0 = ones(240,320);
kdata = kdata.*repmat(sqrt(w),1,1,20);
nft = NUFFTx(k,w,smap_0);
Image = zeros(240,320,20);
for i = 1:20
    Image(:,:,i) = reshape(nft'*kdata(:,:,i),size(smap_0));
end
Image_all = sqrt(sum(conj(Image).*Image,3))+eps;
mask = abs(Image_all)>0.08*max(abs(Image_all(:)));
mask = imdilate(mask,ones(7,7));
smap = Image./repmat(Image_all,1,1,20).*(repmat(mask,1,1,20));
save('smap_dce.mat','smap');