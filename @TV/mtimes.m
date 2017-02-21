function res = mtimes(a,b)

if a.adjoint
    res = tranTV(b);
else
    res = zeros([size(b),2]);
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    bx = b(:,:,:,:,[2:end,end-1]);
    res(:,:,:,:,:,1) = bx - b;
    res(:,:,:,:,:,2) = b;
end

function y = tranTV(b)
y = sum(b,6);
