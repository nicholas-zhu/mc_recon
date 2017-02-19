function [SO_error, SO_grad]=penalties_smoothness_2(O,sizeI,Scaling)
% connection with each other
% O : x, y, dim, t
global penalties;
calgrad = true; 

% normalized 
msizeI=sqrt(sum(sizeI.^2));

% Make penalty independent of image size
O(:,:,1,:)=O(:,:,1,:)/(msizeI/Scaling(1));
O(:,:,2,:)=O(:,:,2,:)/(msizeI/Scaling(2));

% Calculate the thin sheet of metal energy
[PEx,PEgx]=calculate_energy_2d(O(:,:,1,:),calgrad,penalties.Verror2d,penalties.Vgrad2d);
[PEy,PEgy]=calculate_energy_2d(O(:,:,2,:),calgrad,penalties.Verror2d,penalties.Vgrad2d);


function [PE,PEg]=calculate_energy_2d(O,calgrad,Verror2d,Vgrad2d)
sizeI = [size(O,1),size(O,2)];
PE=zeros(sizeI-3);
PEg=zeros(sizeI);

for k = 1:size(O,4)
    for i=1:size(O,1)-3;
        for j=1:size(O,2)-3;
           % Get the control points of one grid cell
           P=O(i+(0:3),j+(0:3),1,k)'; P=P(:)';

           % Calculate the 2D bending energy of thin sheet of metal
           PE(i,j)=sum(sum((P'*P).*Verror2d));

           if(calgrad)
               DP=(Vgrad2d*P(:))';
               % Add penalties from this cell to the total penalty gradients of the control points
               PEg(i+(0:3),j+(0:3),1,k)=PEg(i+(0:3),j+(0:3))+reshape(DP,[4 4])';
           end
        end
    end
end
