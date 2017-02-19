function res = mtimes(a,b)

if a.adjoint
    res = tranLR(b);
else   
    % bx = repmat(mean(b,4),1,1,1,size(b,4));
    temp = movepixels_2d_cpx(b,a.B,2); %cubic interpolation
    % reshape block
    [H,shift,M] = I2H(temp,a.bsize,1);
    % optimization
    H1 = n_norm(H,0.02);
    % back
    temp2 = H2I(H1,a.bsize,shift,M,temp);
    res = movepixels_2d_cpx(temp2,a.F,2);
end

function y = tranLR(b)
y = b;


function [H,shift,M] = I2H(I,bsize,flag_r)
shift = zeros(1,2);
if flag_r
    shift = floor(0.99*rand(1,2).*bsize);
end
Isize = [size(I,1),size(I,2)];
M = floor((Isize-shift)./bsize);
H = zeros(prod(M),prod(bsize),size(I,4));
for i = 1:M(1)
    for j = 1:M(2)
        i_shift = (i-1)*bsize(1)+1+shift(1);
        j_shift = (j-1)*bsize(2)+1+shift(2);
        temp = I(i_shift:i_shift+bsize(1)-1,j_shift:j_shift+bsize(2)-1,:,:);
        H((i-1)*M(2)+j,:,:) = reshape(temp,prod(bsize),size(I,4));
    end
end

function I = H2I(H,bsize,shift,M,I0)

H = reshape(H,[size(H,1),bsize,size(I0,4)]);
I = I0;
for i = 1:M(1)
    for j = 1:M(2)
        i_shift = (i-1)*bsize(1)+1+shift(1);
        j_shift = (j-1)*bsize(2)+1+shift(2);
        I(i_shift:i_shift+bsize(1)-1,j_shift:j_shift+bsize(2)-1,:,:) = H((i-1)*M(2)+j,:,:,:);
    end
end

function H_t = n_norm(H,lambda)
H_t = H;
for i = 1:size(H,1)
    [U, S, V] = svd(squeeze(H(i,:,:)),'econ');
    th = lambda*max(S(:));
    S = SoftThresh(S,th);
    H_t(i,:,:) = U*S*V';
end
