function res = mtimes(a,b)

if a.adjoint
    % ksp -> imag
    b = reshape(b,[size(a.Smap),a.sizeM]);
    res = repmat(a.Mask,1,1,1,a.sizeS,1).*b;
    for dim = 1:3
        res = fftshift(ifft(ifftshift(res,dim),[],dim),dim);
    end
    %res = sum(res,4)*a.nMask;
    res = sum(res.*repmat(conj(a.Smap),1,1,1,1,a.sizeM),4)*a.nMask;
else
    % imag -> ksp
    res = reshape(b,size(a.Mask));
    res = repmat(res,1,1,1,a.sizeS,1).*repmat(a.Smap,1,1,1,1,a.sizeM);
    for dim = 1:3
       res = fftshift(fft(ifftshift(res,dim),[],dim),dim);
    end
    res = res.*repmat(a.Mask,1,1,1,a.sizeS,1)/a.nMask;
end
res = res(:);