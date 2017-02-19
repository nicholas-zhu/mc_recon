function Fxy = movepixels_3d_grid(F1x,B)
% Forward x-y based on x grid, x deformed grid, y target grid
% F: F1->x
% B: Bx->1 - By->1
for i = 1:size(B,4)
    for j = 1:size(B,3)
        Fxy(:,:,:,j,i) =  movepixels_2d_single(B(:,:,:,j,i),F1x(:,:,:,1,i),F1x(:,:,:,2,i),F1x(:,:,:,3,i));
    end
end