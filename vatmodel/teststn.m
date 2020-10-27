function mean_isotropic_conductivity_stn = teststn(vol, cond, side, options)

% Here's how to check the conductivity values of the STN (subthalamic 
% nucleus), expected to be around 0.33, as it is grey matter

try
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'stn.mat']);
catch
    % find STN data
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'atlas_index.mat']);
    stnl = atlases.fv{1, 1}.vertices;
    stnr = atlases.fv{1, 2}.vertices;
    save([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'stn.mat'], 'stnl', 'stnr');
    % load stn points
    load stn
end

% select side
if side == 1
    stn = stnl;
elseif side == 2
    stn = stnr;
end

% find conductivity in the STN
vol_stn = unique(knnsearch(vol.ctr, stn));
cond_stn = cond(vol_stn, :);

% compute anisotropic index and equivalent isotropic conductivity
isotropic_conductivity_stn = zeros([size(cond_stn,1),1]);
aniso_index = zeros([size(cond_stn,1),1]);
for i = 1:size(cond_stn,1)
    isotropic_conductivity_stn(i) = isoCond(cond_stn(i,:));
    S = makeSymmetric(cond_stn(i,:));
    [~,D] = eig(S);
    D = diag(D);
    D = sort(D, 'descend');  % order eigenvalues from the highest to the lowest
    aniso_index(i) = 1 - (D(3) / D(1));
end

% compute summary values for the whole stn volume
mean_cond_stn = mean(cond_stn);
mean_isotropic_conductivity_stn = isoCond(mean_cond_stn);
mean_aniso_index = mean(aniso_index);

% save interesting variables
clear i S D vol_stn
if options.savepath(end) == filesep
    options.savepath(end) = [];
end
save([options.savepath, filesep, 'teststn_output.mat'])

end  % function


function I = isoCond(condArray)
% 
% Returns the equivalent isotropic conductivity I of the array condArray 
% describing the ellipsoid for the anisotropic conductivity
% 

S = makeSymmetric(condArray);
[~,D] = eig(S);
D = diag(D);
I = (prod(D))^(1/3);

end  % subfunction


function S33 = makeSymmetric(array6)

S33 = zeros(3);
S33(1,1) = array6(1);
S33(2,2) = array6(2);
S33(3,3) = array6(3);
S33(1,3) = array6(6);
S33(3,1) = array6(6);
S33(1,2) = array6(4);
S33(2,3) = array6(5);
S33(3,2) = array6(5);
S33(2,1) = array6(4);

end  % subfunction


