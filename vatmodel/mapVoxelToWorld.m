function worldCoordinates = mapVoxelToWorld(voxelCoordinates,image)
voxelCoordinates = [voxelCoordinates';ones(1,size(voxelCoordinates,1))];
worldCoordinates = image.mat*voxelCoordinates;
worldCoordinates = worldCoordinates(1:3,:)';
end