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
        bt = circshift(b,i-1,M_state);
        F = circshift(a.F,i-1,M_state);
        B = circshift(a.B,i-1,M_state);
        B = B-repmat(B(:,:,:,:,1),1,1,1,size(B,M_state));
        Fxy = movepixels_3d_grid(F,B,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = ifftshift(a.w);
        w(1) = 0;
        w_s = sum(w);
        w = repmat(permute(w,[2,3,4,5,1]),size(b,1),size(b,2),size(b,3),1,1);
        bt = movepixels_3d_cpx(bt,Fxy,1);
        bt = sum((bt.*w),4)/w_s;
        bx(:,:,:,i) = bt;
    end
    res(:,:,:,:,:,1) = bx - b;
    res(:,:,:,:,:,2) = b;
end

function y = tranTV(b)
y = sum(b,6);
