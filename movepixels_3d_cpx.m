function Iout = movepixels_3d_cpx(I,B)
for i = 1:size(I,ndims(I))
    Iout(:,:,:,1,i) =  movepixels_3d_single(real(I(:,:,:,i)),B(:,:,:,1,i),B(:,:,:,2,i),B(:,:,:,3,i))+...
        1i*movepixels_2d_double(imag(I(:,:,:,i)),B(:,:,:,1,i),B(:,:,:,2,i),B(:,:,:,3,i));
end