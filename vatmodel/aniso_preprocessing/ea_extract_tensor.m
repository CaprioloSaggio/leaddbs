function dti = ea_extract_tensor(diffusion_image_path, bvals, bvecs)
% 
% INPUT:
% diffusion_image_path: the path to the diffusion image of which we want to
%                       extract the diffusion tensor
% bvals: the path to the text file containing b-value used during the
%        acquisition of dti. It must be 1xN or Nx1 where N is the number of
%        volumes in dti
% bvecs: the path to the text file containing gradient orientations used 
%        during the acquisition of dti. It must be 3xN or Nx3 where N is 
%        the number of volumes in dti
% 
% OUTPUT:
% dti: an LxMxNx6 array representing the diffusion tensor for each voxel.
%      L, M, N are the dimensions of the volumes (N is the number of slices) 
%      and the 4th dimensions contains the elements of the symmetric 
%      positive definite tensor representing diffusivity in each voxel in 
%      the order [xx, yy, zz, xy, yz, xz]
% 
% the core of the code has been downloaded from https://research.dwi.ufl.edu/projects/fanDTasia/tutorial.php
% and then modified.
% 

% load bvecs and check orientation of the array
bvecs = load(bvecs);
if size(bvecs,1)==3
    bvecs=bvecs';
end

% load bvals and check orientation of the array
bvals = load(bvals);
if size(bvals,1)==1
    bvals=bvals';
end

% read the image
% TODO: ask Alba if it makes sense to have negative values in a diffusion image
S = ea_load_nii(diffusion_image_path);
% bias = min(dti,[],'all');
% if bias < 0
%     S = S + bias;
% end
S = abs(S);

% check correspondence between image and acquisition parameters
% TODO: this part will be replaced by a more strict error thrower in a
%       future version
if size(S,4)<size(bvecs,1)
    bvecs = bvecs(1:size(S,4),:);
    warning("bvecs shows a different number of values than the image it is supposed to describe. bvecs will be cut");
elseif size(S,4)>size(bvecs,1)
    error("bvecs has not enough values to describe the diffusion image")
end

if size(S,4)<size(bvals,1)
    bvals = bvals(1:size(S,4));
    warning("bvals shows a different number of values than the image it is supposed to describe. bvals will be cut");
elseif size(S,4 )>size(bvecs,1)
    error("bvals has not enough values to describe the diffusion image")
end

% compute the diffusion tensor with non-negative constraints to get a
% positive definite array for each voxel
G=constructMatrixOfMonomials(bvecs, 2);
C=constructSetOf81Polynomials(2)';
P=G*C;P=[-diag(bvals)*P ones(size(bvecs,1),1)];
dti=zeros(size(S,1),size(S,2),size(S,3), 6); 
for k=1:size(S,3)
    for i=1:size(S,1)
       for j=1:size(S,2)
          y=double(log(squeeze(S(i,j,k,:))));
          x=lsqnonneg(P, y);
          T = C * x(1:81);
          dti(i,j,k,:)=[T(6) T(3) T(1) T(5)/2 T(2)/2 T(4)/2];
       end
    end
end

end