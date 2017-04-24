function res = mtimes(a,b)

if a.adjoint
    res = tranTV(b);
else
    res = zeros([size(b),2]);
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    M_state = 5;
    bx = zeros(size(b));
    for i = 1:size(b,M_state)
        bt = circshift(b,1-i,M_state);
        F = circshift(a.F,1-i,M_state);
        B = circshift(a.B,1-i,M_state);
        Bz = B-repmat(B(:,:,:,:,1),1,1,1,1,size(B,M_state));
        Fxy = movepixels_3d_grid(F,Bz);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = ifftshift(a.w);
        w_s = sum(w);
        w = repmat(permute(w,[2,3,4,5,1]),size(b,1),size(b,2),size(b,3),1,1);
        bt = movepixels_3d_cpx(bt,Fxy);
        bt = sum((bt.*w),5)/w_s;
        bx(:,:,:,1,i) = bt;
    end
    res(:,:,:,:,:,1) = b - bx;
    res(:,:,:,:,:,2) = bx;
end

function y = tranTV(b)
y = sum(b,6);