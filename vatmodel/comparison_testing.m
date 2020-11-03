%% COMPARISON TESTING
tic
%% initialization
% paths
reference_image_path = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\WF\anat_t1.nii';
patient_path = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\WF\stimulations\native';
save_filename = 'unit_testing.txt';
methods = {'horn'; 'anisohorn'};
direction = {'dir', 'ring'};
amplitude = [1 3 5];
unit = {'mA', 'V'};
side = ['r', 'l'];

% structural image and grid
reference_image = ea_load_nii(reference_image_path);
grid.pos = build_grid(reference_image);
grid.tree = KDTreeSearcher(grid.pos);

% result structure
volume_difference = zeros([length(direction), length(amplitude), length(unit), length(side)]);
dice_index = zeros([length(direction), length(amplitude), length(unit), length(side)]);


%% start loop
for i_direction = 1:length(direction)
    for i_amplitude = 1:length(amplitude)
        for i_unit = 1:length(unit)
            for i_side = 1:length(side)
                %% choose stimulation
                stimulation = [direction{i_direction}, '_', num2str(amplitude(i_amplitude)), unit{i_unit}];
                
                if side(i_side)=='r'
                    vat = 'vat_right.mat';
                elseif side(i_side)=='l'
                    vat = 'vat_left.mat';
                else
                    error('Side must be either r (right) or l (left)')
                end
                
                %% load and prepare data
                % horn
                data_horn = load([patient_path,filesep,methods{1},'_',stimulation,filesep,vat]);
                vta_index = unique(knnsearch(grid.tree, data_horn.vatfv.vertices));
                binary_horn = zeros(reference_image.dim);
                binary_horn(vta_index) = 1;

                % aniso
                data_aniso = load([patient_path,filesep,methods{2},'_',stimulation,filesep,vat]);
                vta_index = unique(knnsearch(grid.tree, data_aniso.vatfv.vertices));
                binary_anisohorn = zeros(reference_image.dim);
                binary_anisohorn(vta_index) = 1;
                
                volume_difference(i_side, i_amplitude, i_unit, i_direction) = abs(data_aniso.vatvolume - data_horn.vatvolume) / data_horn.vatvolume;
                dice_index(i_side, i_amplitude, i_unit, i_direction) = dice(binary_horn, binary_anisohorn);
            end
        end
    end
end
toc

% export data to txt file
direction_table = [repmat(direction(1), length(amplitude)*length(unit)*length(side),1); repmat(direction(2), length(amplitude)*length(unit)*length(side),1)];
unit_table = repmat([repmat(unit(1), length(amplitude)*length(side),1); repmat(unit(2), length(amplitude)*length(side),1)], length(direction),1);
amplitude_table = repmat([repmat(amplitude(1), length(side),1); repmat(amplitude(2), length(side),1); repmat(amplitude(3), length(side),1)], length(direction)*length(unit),1);
side_table = repmat([side(1); side(2)], length(direction)*length(unit)*length(amplitude),1);
volume_difference_table = reshape(volume_difference, [], 1);
dice_index_table = reshape(dice_index, [], 1);
T = table(direction_table, unit_table, amplitude_table, side_table, volume_difference_table, dice_index_table);
writetable(T,[patient_path,filesep,save_filename]);
