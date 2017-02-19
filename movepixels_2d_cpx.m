function Iout = movepixels_2d_cpx(I,B,flag)
for i = 1:size(I,ndims(I))
    Iout(:,:,1,i) =  movepixels_2d_double(real(I(:,:,i)),B(:,:,1,i),B(:,:,2,i),flag)+...
        1i*movepixels_2d_double(imag(I(:,:,i)),B(:,:,1,i),B(:,:,2,i),flag);
end