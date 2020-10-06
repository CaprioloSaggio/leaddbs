function [stiff, diinsy, cols, sysmat] = ea_sb_stiff_matrix(cond, vol)

% based on SB_CALC_STIFF found in the FieldTrip package
%
% 
% INPUT
% cond: conductance tensor represented as a Nx6 array of conductivity
%       values
% vol: volume mesh structure containing nodes in the field called 'pos' and
%      the elements connectivity in a field called 'tet'. It must have N
%      elements
% 
% 
% OUTPUT
% stiff: sparse matrix containing in the element (i,j) the conductivity
%        between nodes i and j
% diinsy:
% cols:
% sysmat:

% TODO: check if I need all those outputs

node = double(vol.pos);

npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);  % ##### original
%         mele = size(vol.tet,1)*ones(1,6);
        elem = vol.tet;  % ##### original
%         elem = repmat(vol.tet, 1, 6);
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
%         mele = size(vol.tet,1)*ones(1,6);
        elem = vol.tet';
%         elem = repmat(vol.tet', 1, 6);
    else
        error('The volume structure in input has tet (tetrahedrals) field with wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
else
    error("Structure representing the tetrahedral mesh must have a 'tet' field where list of nodes belonging to an element is stored")
end

if size(vol.tet, 1) == size(cond, 1)
    cond = double(cond);
else
    error("The number of conductivity values doesn't match the number of elements in the mesh")
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider swapping node 3 and 4');
    end
end

try
    [diinsy,cols,sysmat] = ea_calc_stiff_matrix_val_wrapper(node,elem,cond,mele);  % #####
    [diinsy,cols,sysmat] = calc_stiff_matrix_val(node,elem,cond,mele);  % ##### original
    ea_delete([pwd, filesep, 'fort.6']);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));
end
