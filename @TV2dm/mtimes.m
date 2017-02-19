function res = mtimes(a,b)

if a.adjoint
    res = tranTV(b);
else
    res = zeros([size(b),2]);
    % could design weights for each bin
    % bx = repmat(mean(b,4),1,1,1,size(b,4));
    bx = zeros(size(b));
    for i = 1:size(b,4)
        bt = circshift(b,i-1,4);
        F = circshift(a.F,i-1,4);
        B = circshift(a.B,i-1,4);
        B = B-repmat(B(:,:,:,1),1,1,1,size(B,4));
        Fxy = movepixels_2d_grid(F,B,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = ifftshift(a.w);
        w(1) = 0;
        w_s = sum(w);
        w = repmat(permute(w,[2,3,4,1]),size(b,1),size(b,2),1,1);
        bt = movepixels_2d_cpx(bt,Fxy,1);
        bt = sum((bt.*w),4)/w_s;
        bx(:,:,:,i) = bt;
    end
    res(:,:,:,:,1) = b - bx;
    res(:,:,:,:,2) = bx;
end

function y = tranTV(b)
y = sum(b,5);
