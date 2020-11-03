function [mean_snr_vol, mean_snr_stn] = testsnr(vol, mri, side, options)
% 
% INPUT
% vol: vol structure of ea_genvat_aniso, which contains pos, tet, ctr...
% mri: structural image read by using ea_load_nii
% side: side of stimulation (1,right or 2,left)
% options: the options structure in any ea_genvat_*
% 
% OUTPUT
% mean_snr_vol: mean SNR in the whole region taken as subset of the head
%               for field computation
% mean_snr_stn: mean SNR in the STN
% 
% 

% build grid
grid = build_grid(mri);
% build the KDtree
gridtree = KDTreeSearcher(grid);

snr = niftiread([options.root, options.patientname, filesep, 'snr_filtered.nii.gz']);

% find mean SNR in the whole volume
vol_index = knnsearch(gridtree, vol.ctr);
snr_vol = snr(vol_index);
mean_snr_vol = mean(snr_vol);

% find mean SNR in the STN only
try
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'stn.mat']);
catch
    % find STN data
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'atlas_index.mat']);
    stnl = atlases.fv{1, 1}.vertices;
    stnr = atlases.fv{1, 2}.vertices;
    save([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'stn.mat'], 'stnl', 'stnr');
    % load stn points
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'stn.mat']);
end
% select side
if side == 1
    stn = stnl;
elseif side == 2
    stn = stnr;
end
stn_index = knnsearch(gridtree, stn);
snr_stn = snr(stn_index);
mean_snr_stn = mean(snr_stn);






