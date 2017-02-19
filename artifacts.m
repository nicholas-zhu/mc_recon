% artifacts value evaluation
clc;
clear;
addpath(genpath('ESPIRiT'));
addpath(genpath('fminlbfgs_version2c'));
addpath ~/bart/matlab/
load nocart_data.mat;
load smap.mat

% binning, how many bins
Res_d = diff(filter([0.1 0.2 0.4 0.2 0.1],1,Res_Signal));
Res = sign([Res_d(1);Res_d]).*Res_Signal;
motion_signal = [Res,(1:length(Res_Signal))'];
motion_signal = sortrows(motion_signal,1);
motion_signal = motion_signal(:,2);
nbin = 10;

k_sort = double(k(:,motion_signal));
kdata_sort = double(kdata(:,motion_signal,:));
b1 = smap;
w = sqrt(abs(k_sort));
k_full = k_sort;
E0 = NUFFTx(k_sort,w,double(b1));
kdata_full = kdata_sort.*repmat(w,1,1,size(b1,3),1);
I0 = reshape(E0'*kdata_full,240,320);

k_sort = reshape(k_sort,size(k,1),size(k,2)/nbin,1,nbin);
kdata_sort = permute(reshape(kdata_sort,size(k,1),...
    size(k,2)/nbin,nbin,size(kdata,3)),[1 2 4 3]);
N1 = 1;
N2 = 1;
w = sqrt(abs(k_sort));
k_sort = real(k_sort)/N1+1i*imag(k_sort)/N2;
E1 = NUFFTx(k_sort,w,double(b1));
kdata_sort = kdata_sort.*repmat(w,1,1,size(b1,3),1);
I1 = reshape(E1'*kdata_sort,240,320,1,nbin);