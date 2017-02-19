% real value intensity based dynamic image registration work flow

% dicom file load
N = 14;
for i = 1:N
    I(:,:,i) = dicomread(['E4288S14l',num2str(i),'.DCM']);
end
mask = (mean(I,3)>(0.4*mean(I(:))));
se = ones(3);
mask = imerode(mask,se);
mask = imdilate(mask,se);
%% normalization & intensity adjustment
gamma = 1;
I = single(I);
I = I./max(I(:));
I = reshape(imadjust(I(:,:),[0,1],[0,1],gamma),size(I));

% affine registration??
Options.Similarity = 'sd';
Options.Spacing = [32 32];
Options.Penalty = 1e-2;
Options.Registration = 'NonRigid';

% for i = 1:14
%     [Ireg(:,:,i),O_trans,Spacing,M,B(:,:,:,i),F(:,:,:,i)] = image_registration(I(:,:,i),I(:,:,1),Options);
% end

%% continuous registration
% 
spacing = [32 32];
IterM = min(floor(log2(spacing/4))+1);
center = size(I(:,:,1))/2;
O_trans = make_init_grid_2(spacing,size(I(:,:,1)),[],N);

iter = 1;
Hsize=round(0.25*(size(I,1)/size(O_trans,1)+size(I,2)/size(O_trans,2)));
for i = 1:N
    I_1(:,:,i) = imfilter(I(:,:,i),fspecial('gaussian',[Hsize Hsize],Hsize/5));
end

resize_per = 2^(iter-1);
MASK = imresize(mask,1/resize_per);
I_r = zeros([size(MASK),N]);
for i = 1:N
    I_r(:,:,i) =imresize(I_1(:,:,i),1/resize_per);
end

Spacing_small=spacing/resize_per;

if(size(I1,3)<4)
    for i = 1:N
        [O_error(:,:,:,i),O_grad(:,:,:,i)] = bspline_error_2d_double(double(O_grid(:,:,1,i)),...
            double(O_grid(:,:,2,i)),double(I_1(:,:,:,i)),double(I_1(:,:,:,1)),double(Spacing_small(1)),double(Spacing_small(2)));
    end
    [SO_error,SO_grad] = penalties_smoothness_2(O_grid,size(I1),options.scaling);
end

disp(['Error' num2str(O_error)]);
disp(['Smoothness Error' num2str(SO_error)]);
disp(['Abs Mean Error Gradient' num2str(mean(abs(O_grad(:))))]);
disp(['Abs Max  Error Gradient' num2str(max(abs(O_grad(:))))]);
disp(['Abs Max  Smoothness Gradient' num2str(max(abs(SO_grad(:))))]);
disp(['Abs Mean Smoothness Gradient' num2str(mean(abs(SO_grad(:))))]);

O_grad=O_grad+SO_grad*penaltypercentage;