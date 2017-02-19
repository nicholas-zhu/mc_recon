function y = mft(img,M,S)
img = reshape(img,size(M));
res = repmat(img,1,1,1,size(S,4),1).*repmat(S,1,1,1,1,size(img,5));
for dim = 1:3
    y = fftshift(ifft(ifftshift(res,dim),[],dim),dim);
end
y = y.*repmat(M,1,1,1,size(S,4),1);
y = y(:);