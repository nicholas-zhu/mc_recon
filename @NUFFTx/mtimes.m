function ress = mtimes(NFT,b)

% motion nufft
if NFT.adjoint
    b = reshape(b,NFT.dataSize);
    ress = zeros(NFT.imSize2);
    res = zeros(size(NFT.S));
    for i = 1:size(b,4)
        for j = 1: size(b,3)
            bx = b(:,:,j,i).*NFT.w(:,:,:,i);
            res(:,:,j) = reshape(nufft_adj(bx(:),NFT.st{i})/sqrt(prod(NFT.imSize)),NFT.imSize(1),NFT.imSize(2));
        end
        ress(:,:,1,i)=sum(res.*conj(NFT.S),3)./(sum(abs((squeeze(NFT.S))).^2,3)+eps);
    end
    ress = ress(:)*size(NFT.w,1)*pi/2/size(NFT.w,2);
else
    b = reshape(b,NFT.imSize2);
    ress = zeros(NFT.dataSize);
    for i=1:size(b,4)
        for j=1:size(NFT.S,3)
            bx=b(:,:,1,i).*NFT.S(:,:,j); 
            ress(:,:,j,i) = reshape(nufft(bx,NFT.st{i})/sqrt(prod(NFT.dataSize(1:2))),NFT.dataSize(1),NFT.dataSize(2)).*NFT.w(:,:,i);
        end
    end
    ress = ress(:);
end