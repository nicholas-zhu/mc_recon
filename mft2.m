function y_i = mft2(img,M,S)
img = reshape(img,size(M));
res = repmat(img,1,1,1,size(S,4),1).*repmat(S,1,1,1,1,size(img,5));
for dim = 1:3
    y = fftshift(ifft(ifftshift(res,dim),[],dim),dim);
end
y = y.*repmat(M,1,1,1,size(S,4),1);

for dim = 1:3
    y = fftshift(fft(ifftshift(y,dim),[],dim),dim);
end

y_i = sum(y.*repmat(conj(S),1,1,1,1,size(M,5)),4);
y_i = y_i(:);