function res = mtimes(a,b)

if a.adjoint
    res = tranTV(a.F,b);
else
    res = zeros([size(b),2]);
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    b = movepixels_2d_cpx(b,a.B,1);
    bx = b(:,:,:,[2:end,end-1]);
    res(:,:,:,:,1) = b - bx;
    res(:,:,:,:,2) = bx;
end



function y = tranTV(F,b)
y = sum(b,5);
for i =1:size(y,4)
    y = movepixels_2d_cpx(y,F,1);
end