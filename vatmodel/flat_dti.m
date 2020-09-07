function dti_flat = flat_dti(dti)

% flattens a dti 4D volume in a 2D array Nx6, where N is the number of the
% voxels in the DTI.
volume_voxels = prod(size(dti,1:3));
dti_flat = zeros(volume_voxels, 6);

for k = 1:size(dti,3)
    delta_slice = (k-1)*size(dti,2)*size(dti,1);
    for i=1:size(dti,1)
        delta_row = (i-1)*size(dti,2);
        delta = delta_slice + delta_row;
        
        dti_flat(delta+1:delta+size(dti,2),:) = dti(i,:,k,:);
    end
end

end  % function

% volume_voxels = prod(size(dti,1:3));
% 
% dti_flat = zeros(volume_voxels, 6);
% 
% for i=size(dti,1):-1:1
%     delta_row = (i-1)*size(dti,3)*size(dti,2);
%     
%     for k=size(dti,3):-1:1
%         delta_slice = (k-1)*size(dti,2);
%         delta = volume_voxels - delta_slice - delta_row;
%         
%         dti_flat(delta:-1:1+delta-size(dti,2),:) = dti(i,:,k,:);
%     end
% end
% 
% end  % function

