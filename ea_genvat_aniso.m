function varargout = ea_genvat_aniso(varargin)
% 
% this function is called at line 1059 of function ea_stimparams.m in the
% following fashion:
% [stimparams(1,side).VAT(el).VAT,volume]=feval(ea_genvat,elstruct(el).coords_mm,getappdata(handles.stimfig,'S'),side,options,stimname,options.prefs.machine.vatsettings.aniso_ethresh,handles.stimfig);
%
% TODO: clean the code from commented lines
% TODO: make the load of the dti tensor clean (now a variable called 'DTI' is loaded and then substituted by one called 'dti')
% TODO: make it monopolar, or make it work with more than one contact
%       stimulating
% TODO: make the substitution of the electrode conductivity faster
% TODO: test w/ and w/o block at line 664
% 


if nargin==5
    coords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    thresh=options.prefs.machine.vatsettings.aniso_ethresh;  % standard 0.2 V/m;
    
elseif nargin==7
    coords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    thresh=varargin{6};
    lgfigure=varargin{7};

elseif nargin==1  % in this case I just return the name of the method
    if ischar(varargin{1}) 
        varargout{1} = 'Anisotropic';
        return
    end
end

% check if the VTA must be estimated for the side in scope
if not(any(S.amplitude{1,side}))
    varargout{1}=nan;
    varargout{2}=nan;
    varargout{3}=nan;
    fprintf('SIDE %d HAS NOT BEEN STIMULATED \n', side);
    return
end

clear varargin

%% debug settings
dbg_vis = 0;  % flag to state if data visualization for debugging is to be run
dbg_recompute = 0;  % flag to state if the model has to be recomputed all the times (1) or not (0)
dbg_fast = 0;  % flag which if set to 0 makes the model faster to be computed, but less accurate
dbg = 0;  % enter in debug mode if set to 1

if dbg; tic; end
recompute = options.overwriteapproved || dbg_recompute;  % the user can also use the checkbox in the LEAD-DBS menu to recompute the model all the times


%% conventions and settings
encaps = options.prefs.machine.vatsettings.aniso_encaps;  % flag to determine if we add the encapsulation (scar) tissue to the model
SI = 0;  % flag to determine if to use the International System or not, i.e. if to use m (1) or mm (0). Whenever this flag is changed, the model must be recomputed
tensor_name = get(options.prefs.machine.vatsettings.aniso_tensor_name,'String');
dti_tensor = tensor_name;  % name of the file containing the diffusion tensor, "dti_tensor.mat default

% check if the diffusion tensor exists in patient folder and has a supported data structure (.mat ot .nii)
if not(exist([options.root,options.patientname,filesep,dti_tensor], 'file') == 2)
    error("File %s not found in the patient folder, make sure you have named it properly", dti_tensor)
end
if not(isequal(dti_tensor(end-3:end),'.nii') || isequal(dti_tensor(end-3:end),'.mat'))
    error("File containing the diffusion tensor must have .nii or .mat file format. Make sure you specified it properly in Settings");  % check the user inserted the mat
end

if SI
    thresh = thresh * 1e3;
end

%% ________________________________________________________________________
%% ELECTRODE DATA
ea_dispt('Building electrode model...')


%% get the electrode model and the relative patient-specific transformation
% Important to load in reco from a new since we need to decide whether to
% use native or template coordinates. Even when running in template space,
% the native coordinates are sometimes used (VTA is then calculated in 
% native space and ported to template).

resultfig=getappdata(lgfigure,'resultfig');
options.loadrecoforviz=1;
[coords_mm,trajectory,markers]=ea_load_reconstruction(options);
elstruct(1).coords_mm=coords_mm;
elstruct(1).coords_mm=ea_resolvecoords(markers,options);
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;

elspec=getappdata(resultfig,'elspec');
coords=coords{side};
setappdata(resultfig,'elstruct',elstruct);

% ##########
% coords_m{1, 1} = coords_mm{1, 1} / 1e3;  % from mm to m
% coords_m{1, 2} = coords_mm{1, 2} / 1e3;

% Add stretchfactor to elstruct simply for purpose of checking if headmodel
% changed. Larger stim amplitudes need larger bounding boxes so
% stretchfactor must be incorporated here.
if max(S.amplitude{side})>4
    elstruct.stretchfactor=0.75;
else
    elstruct.stretchfactor=0.5;
end

% load and transform volumetric mesh of the electrode (w/o information about 
% insulation or contact)
elmodel_path=[ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'];
elmodel = load(elmodel_path);

if dbg_vis
    figure
    plot3(elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3), 'b.');
end


%% ________________________________________________________________________
%% BUILD MESH
    
% if it already exists, load the mesh
    if side==1
        rl = 'right';  % rl (standing for right/left) is a flag that explicit the side considered at the present call to the function
    elseif side==2
        rl = 'left';
    else
        error("The considered side is neither 1 (right) or 2 (left). Check the call to ea_genvat_aniso function is done properly")
    end
if (exist([options.root, options.patientname, filesep, 'stimulations', filesep, 'mesh_', rl, '.mat'], 'file') == 2) && not(recompute)
    if dbg
        load([options.root, options.patientname, filesep, 'stimulations', filesep, 'mesh_', rl, '.mat'], 'cylinder', 'T', 'vol', 'electrode', 'contacts_vertices', 'encaps_index', 'elmodel', 'elrad')
    else
        load([options.root, options.patientname, filesep, 'stimulations', filesep, 'mesh_', rl, '.mat'], 'vol', 'T', 'electrode', 'contacts_vertices', 'encaps_index', 'elmodel', 'elrad')
    end
% ortherwise compute it 
else
    %% ________________________________________________________________________
    %% BOUNDING CYLINDER
    fprintf('\nVAT ESTIMATION ON SIDE %i \n', side)
    ea_dispt('Building the bounding cylinder...')

    if max(S.amplitude{side})>4
        cylinder.stretchfactor=0.75*(max(S.amplitude{side})/2.5);
    else
        cylinder.stretchfactor=0.5;
    end

    % define a bounding cylinder to which restrict the mesh
    cylinder.upper_bound = 15+(20*cylinder.stretchfactor);
    cylinder.lower_bound = -20*cylinder.stretchfactor;
    cylinder.vol_top = [0, 0, cylinder.upper_bound];
    cylinder.vol_bottom = [0, 0, cylinder.lower_bound];
    cylinder.cyltrisize = 10;                 % maximum triangle size of the bounding cyl
    cylinder.cyltetvol = 10;                  % maximum tetrahedral volume in the mesh
    cylinder.cylradius = 40*cylinder.stretchfactor;   % define the radius of the bounding cylinder
    cylinder.ndiv=50;                        % division of circle for the bounding cylinder

    % define and transform the cylinder, directly obtaining the mesh
    [vol.pos,vol.face,vol.tet] = meshacylinder(cylinder.vol_bottom, cylinder.vol_top, cylinder.cylradius, cylinder.cyltrisize, cylinder.cyltetvol, cylinder.ndiv);
    vol.tet = vol.tet(:,1:4);  % eliminate last column, that is the one used for homogeneous representation
    
    
    %% ________________________________________________________________________
    %% mesh refinement with nodes of the electrode    
    ISO2MESH_SURFBOOLEAN='cork';
%     for attempt=1:4
%         try 
            [vol.pos, vol.face] = surfboolean(vol.pos, vol.face(:,1:3), 'resolve', elmodel.node, elmodel.face);  % add electrode mesh
            [vol.pos, vol.face] = meshcheckrepair(vol.pos, vol.face, 'dup');  % clean the mesh
            [vol.pos, vol.face] = meshcheckrepair(vol.pos, vol.face, 'deep');
            [vol.pos, vol.face] = meshcheckrepair(vol.pos, vol.face, 'dup');
            [vol.pos, vol.tet, vol.face] = s2m(vol.pos, vol.face(:,1:3), 1, 1);  % obtain the definitive mesh
%             break
%         catch  % try to move the points
%             fprintf("\nAttempt %i failed in generating the mesh\n", attempt);
%             vol.pos = vol.pos + randn/700;
%             elmodel.node = elmodel.node + randn/700;
%             if attempt==4
%                 error("Not able to generate the mesh, check the data of the electrode")
%             end
%         end
%     end
    clear ISO2MESH_SURFBOOLEAN
    vol.tet = vol.tet(:,1:4);
    vol.ctr = meshcentroid(vol.pos, vol.tet);
    
    
    %% ________________________________________________________________________
    %% BOUNDARY NODES
    % find boundary points in the volume of interest
    vol.boundary = unique(boundary(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3)));
    if dbg_vis
        figure 
        plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'b.')
        hold on
        plot3(vol.pos(vol.boundary, 1), vol.pos(vol.boundary, 2), vol.pos(vol.boundary, 3), 'r.')
        title('boundary nodes')
        legend('nodes inside the volume of interest','boundary nodes', 'location', 'northeast')
    end
    

    %% ________________________________________________________________________
    %% electrode transformation and contacts
    % compute transformation from general to patient specific electrode model
    % (surface, containing info on insulation or contact)
    [~,~,T,electrode]=ea_buildelfv(elspec,elstruct,side);  % T is the transformation between model space and patient space

    
% ##########
%     %% ________________________________________________________________________
%     %% cut the electrode at the size of the bounding cylinder
%     if not(refine)
%         elmodel.node(elmodel.node(:,3)>cylinder.upper_bound | elmodel.node(:,3)<cylinder.lower_bound, :) = [];
% 
%         if dbg_vis
%         %     elmodel.ctr = tetctr(elmodel.node, elmodel.face+1);  % find the centroids of the electrode elements
%             figure
%             plot3(elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3), 'bx');
%         %     plot3(elmodel.ctr(:,1), elmodel.ctr(:,2), elmodel.ctr(:,3), 'bx');
%             hold on
%             plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'ro')
%         end
%     end


    %% ________________________________________________________________________
    %% ENCAPSULATION LAYER
    % it is modelled assuming it as a 0.5mm thick layer around the whole length of
    % the electrode lead
    vol.r = vecnorm(vol.ctr(:,1:2)')';  % find radial distance of each point in the bounding cylinder from the center of the electrode
    elrad = elspec.lead_diameter / 2;  % find radius of the electrode lead
    encapsulation_thickness = 0.5;  % 0.5 mm according to Gunalan et al 2017
    encaps_index = find(vol.r>elrad & vol.r<(elrad+encapsulation_thickness) & vol.ctr(:,3)>min(elmodel.node(:,3)));  % find all elements in the cylinder mesh that correspond to encapsulation tissue
    encaps_index = [encaps_index; find(vol.r<elrad & vol.ctr(:,3)>min(elmodel.node(:,3))-encapsulation_thickness & vol.ctr(:,3)<min(elmodel.node(:,3)))];  % rough modelling for the curved contact 0
        
    if dbg_vis
        figure
        plot3(vol.ctr(:,1), vol.ctr(:,2), vol.ctr(:,3), 'r.')
        hold on
        plot3(vol.ctr(encaps_index,1), vol.ctr(encaps_index,2), vol.ctr(encaps_index,3), 'b.')
        legend('cylindrical volume nodes', 'fibrotic tissue nodes', 'location', 'northeast')
        title('encapsulation layer')
    end

    
        %% ________________________________________________________________________
    %% transform the electrode model into patient space
    elmodel.node = T*[elmodel.node, ones(size(elmodel.node,1),1)]';
    elmodel.node = elmodel.node(1:3,:)';

% ##########
%     % put together the mesh containing all the contacts in the electrode
%     contacts_vertices = [];
%     for i=1:length(electrode.contacts)
%         contacts_vertices = [contacts_vertices; electrode.contacts(i).vertices]; 
%     end


    %% ________________________________________________________________________
    %% BOUNDING CYLINDER (transform)
    % transform using the same transformation of the electrodes (from model to patient space)
    vol.pos = T*[vol.pos, ones(length(vol.pos), 1)]';
    vol.pos = vol.pos(1:3,:)';
    vol.ctr = T*[vol.ctr, ones(length(vol.ctr), 1)]';
    vol.ctr = vol.ctr(1:3,:)';
    vol.unit = 'mm';

    if dbg_vis
        figure
        title('Bounding cylinder');
        plot3(vol.pos(:,1),vol.pos(:,2),vol.pos(:,3),'r*');
        hold on
        plot3(contacts_vertices(:,1), contacts_vertices(:,2), contacts_vertices(:,3), 'ro')
        plot3(elmodel.node(:,1),elmodel.node(:,2),elmodel.node(:,3),'g.');
    end

    
    %% find boundary points in the volume of interest
    vol.boundary = unique(boundary(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3)));
    
    if dbg_vis
        figure 
        plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'b.')
        hold on
        plot3(vol.pos(vol.boundary, 1), vol.pos(vol.boundary, 2), vol.pos(vol.boundary, 3), 'r.')
        title('boundary nodes')
        legend('nodes inside the volume of interest','boundary nodes', 'location', 'northeast')
    end

    
    
    %% ________________________________________________________________________
    %% SAVE VOLUME AND RELATED VARIABLES    
    save([options.root, options.patientname, filesep, 'stimulations', filesep, 'mesh_', rl, '.mat'], 'vol', 'cylinder', 'T', 'electrode', 'encaps_index', 'elmodel', 'elrad')

end

    
%% ________________________________________________________________________
%% HEADMODEL
ea_dispt('Starting FEM headmodel generation...')

if side==1
    rl = 'right';
elseif side==2
    rl = 'left';
else
    error("The considered side is neither 1 (right) or 2 (left). Check the call to ea_genvat_aniso function is done properly")
end
if (exist([options.root, options.patientname, filesep, 'stimulations', filesep, 'conduction_model_', rl, '.mat'], 'file') == 2) && not(recompute)
    load([options.root, options.patientname, filesep, 'stimulations', filesep, 'conduction_model_', rl, '.mat'], 'vol', 'cond', 'el_cond')
else
    %% read image
    ea_dispt('Reading images...')
    
    anat = dir([options.root, options.patientname, filesep, options.prefs.prenii_searchstring]);
    anat = [options.root,options.patientname,filesep,anat(1).name];  % anat is the path to the file called (usually) anat_t1.nii in the patient folder

    % load structural image
    mri = ea_load_nii(anat);
    
    % load diffusion tensor
    if isequal(dti_tensor(end-3:end),'.nii')
        dti = niftiread(dti_tensor);
    else
        load([options.root,options.patientname,filesep,'dti_tensor.mat'],'dti');
    end
    

    %% build the grid of points corresponding to the voxels of the volume
    % check if the dti has the same dimensions of the anatomical image
    if not(isequal(size(dti, 1:3), mri.dim))
        error("The diffusion tensor dti_tensor.nii needs to have the same dimensions of the anatomical image used as anchor image. Perform image warping")
    end

    % build grid from volume 
    mesh_grid.pos = build_grid(mri);
    
    if not(dbg)
        clear mri anat
    end

    
    %% find correspondence with conductivity values
    ea_dispt('Defining the conductivity tensor...')

    % build the KDtree
    mesh_grid.treepos = KDTreeSearcher(mesh_grid.pos);

    % run Nearest Neighbours to get correspondence between mesh elements 
    % and image voxels
    cond_index = knnsearch(mesh_grid.treepos, vol.ctr); % contains the index of the 
                                                        % conductivities that have
                                                        % a match in the mesh                                                
    if not(dbg)
        clear mesh_grid
    end

    
    %% build headmodel
    % get the diffusion data and convert them to conduction data
    raw_cond_tensor = dti * 0.736;  % S*s/mm^2, McIntyre et al 2004
%     raw_cond_tensor = (dti - 0.124e-6) * 0.844;  % S*s/mm^2, Tuch et al 2001 (this formula doesn't guarantee for positive definiteness of the resulting conductivity tensor)
    clear dti

    % rearrange conductivity tensor
    cond_tensor = zeros([prod(size(raw_cond_tensor, 1:3)), 6]);
    for i=1:6
        cond3d = raw_cond_tensor(:,:,:,i);   % this passage is added for more readibility
        cond_tensor(:, i) = reshape(cond3d, [numel(cond3d), 1]);
    end
    
    % select only the conductivities that match elements in the tetrahedral mesh
    cond = cond_tensor(cond_index, :);

% ##########
    % rearrange columns of the conductivity as simbio wants them in the order
    % xx, yy, zz, xy, yz, xz while fsl gives in output xx xy xz yy yz zz
%     cond = [cond(:,1), cond(:,4), cond(:,6), cond(:,2), cond(:,5), cond(:,3)];  #####
%     cond(:,1:3) = abs(cond(:,1:3));  #####
     
    if not(dbg)
        clear cond3d cond_tensor cond_index raw_cond_tensor i
    end
    
    if dbg
        isotropic_conductivity_stn = teststn(vol, cond, side, options);  % it is returned in S/mm
        fprintf("\nDEBUG ________________________________________________________\n")
        fprintf("The average isotropic conductivity in the STN is %f S/m\n", isotropic_conductivity_stn*1e3)
        fprintf("This value is only reliable in native space\n");
        fprintf("______________________________________________________________\n\n\n")
    end


    %% find correspondence with electrodes and include them into the model
    % I set the conductivity of all voxels being close to the centroid of an
    % electrode element to the insulating value of 10^-16 S/m and then I apply the
    % contact conductivity of 10^8 S/m to all the voxels near to the centroid
    % of a contact element. Values are taken from Horn et al 2017
    cond_contacts = options.prefs.machine.vatsettings.aniso_cond_contacts;  % 1e8 * 1e-3 S/mm default
    cond_insulation = options.prefs.machine.vatsettings.aniso_cond_insulation;  % 1e-16 * 1e-3 S/mm default
    cond_encapsulation = options.prefs.machine.vatsettings.aniso_cond_encapsulation;  % according to Gunalan et al 2017 it can be 0.05±0.2 S/m

    % el_cond = unique(knnsearch(vol.ctr, elmodel.ctr));
    if dbg_fast
        el_cond = unique(knnsearch(vol.ctr, elmodel.node, 'K', 5));  % index of all the elements in cylinder corresponding to electrode contacts   
    else
% ##########        
%         el_cond = [];
%         for i=1:length(electrode.insulation)
%             el_cond = [el_cond; find(insurface(electrode.insulation(i).vertices*1e3, electrode.insulation(i).faces, vol.ctr*1e3))];
%         end
        el_cond = find(insurface(elmodel.node, elmodel.face(:,1:3), vol.ctr));
        el_cond = unique(el_cond);
    end
    cond(el_cond,:) = repmat([repmat(cond_insulation,1,3) zeros(1,3)], length(el_cond), 1);  %#ok<FNDSB> % insulation conductivity (isotropic)

    % add encapsulation layer (fibrotic tissue forming around the electrodes)
    if encaps
        cond(encaps_index,:) = repmat([repmat(cond_encapsulation,1,3) zeros(1,3)], length(encaps_index),1);
    end
    
    % find contact elements
    if dbg_fast   
        con_cond = unique(knnsearch(vol.ctr, contacts_vertices));
    else
        inside_active_contacts = false(size(vol.ctr,1),1);
        for c=1:length(electrode.contacts)
            inside_active_contacts = [inside_active_contacts | insurface(electrode.contacts(c).vertices*1e3, electrode.contacts(c).faces, vol.ctr*1e3)]; % #####
        end
% ##########     
%         active_contacts = find(S.activecontacts{1,side});
%         contact_source = zeros(length(active_contacts),1);
%         for c = active_contacts  % 1:length(electrode.contacts)
%             inside_active_contacts = [inside_active_contacts | insurface(electrode.contacts(c).vertices*1e3, electrode.contacts(c).faces, vol.ctr*1e3)]; % #####
% %             inside_active_contacts = inside_active_contacts | inpolyhedron(electrode.contacts(c), vol.ctr); % #####
% 
% %             active_contacts_coords = electrode.contacts(c).vertices;  % in mm  % #####
% %             vol.active = unique(knnsearch(vol.pos, active_contacts_coords));  % #####
% %             contact_source = ea_find_elec_center(vol.active(c), vol.pos);  % #####
% %             inside_active_contacts(inside_active_contacts==contact_source) = [];  % #####
% %             clear vol.active  % #####
% 
% %             active_pos = unique(knnsearch(vol.pos, electrode.contacts(c).vertices));
% %             center_id = ea_find_elec_center(active_pos, vol.pos);
% %             inside_active_contacts(any((vol.tet==center_id)')) = 0;
% 
%         end
        con_cond = find(inside_active_contacts);
    end
    cond(con_cond,:) = repmat([repmat(cond_contacts,1,3) zeros(1,3)], length(con_cond), 1);  % contact conductivity (isotropic)  
    
    if SI
        % convert from mm to m
        cond = cond * 1e3;
        vol.pos = vol.pos * 1e-3;
        vol.ctr = vol.ctr * 1e-3;
        vol.r = vol.r * 1e-3;
        vol.unit = 'm';
    end

    if dbg_vis
        if SI; scale = 1e3; else; scale = 1; end
        figure; 
        plot3(vol.ctr(el_cond,1)*scale, vol.ctr(el_cond,2)*scale, vol.ctr(el_cond,3)*scale, 'b.')
        hold on
        plot3(vol.ctr(con_cond,1)*scale, vol.ctr(con_cond,2)*scale, vol.ctr(con_cond,3)*scale, 'g.')
        plot3(vol.ctr(encaps_index,1)*scale, vol.ctr(encaps_index,2)*scale, vol.ctr(encaps_index,3)*scale, 'r.')
        title('non-DTI nodes')
        legend('insulation', 'contact', 'scar tissue', 'location', 'northwest')
    end

    
    %% ________________________________________________________________________
    %% COMPUTE CONDUCTION MODEL
    ea_dispt('Computing the conduction model...')

% ##########
%     test = cond;
%     cond = cond(:,1);  % #####
    
    try
        vol.stiff = ea_sb_stiff_matrix(cond, vol);
% ##########
%         vol.stiff = simbio_stiff_matrix(cond, vol);
    %     vol.stiff = sb_calc_stiff(cond, vol);
    catch
        vol.tet(:, [3, 4]) = vol.tet(:, [4, 3]);  % necessary not to get 
                                                  % an error from sb_calc_stiff 
                                                  % relative to orientation
        vol.stiff = ea_sb_stiff_matrix(cond, vol);
% ##########
%         vol.stiff = simbio_stiff_matrix(cond, vol);
    %     vol.stiff = sb_calc_stiff(cond, vol);
    end
    vol.method = 'simbio'; 

    % check if the stiffness matrix has all the elements in the diagonal
    % positive, and in case, substitute null values with a very small one for
    % the sake of solvability of the Cholesky factorization during potential
    % computation. Similar approach with negative values on the diagonal, that
    % are considered in their absolute value.
    stiff_diag = full(diag(vol.stiff));
    null_diag_stiff = find(stiff_diag==0);
    if not(isempty(null_diag_stiff))
        for i=1:length(null_diag_stiff)                          % I had to arrange the check this way, 
            j = null_diag_stiff(i);                              % because sometimes MATLAB gives  
            vol.stiff(j,j) = eps;                                % problems of memory outage with 8GB 
        end                                                      % RAM
    end
    
    negative_diag_stiff = find(stiff_diag<0);
    if not(isempty(negative_diag_stiff)) && not(dbg)
        warning("Error in the conductivity values, some of the conductivity values are negative")
    end
    if not(isempty(negative_diag_stiff))
        for i=1:length(negative_diag_stiff)                         % I had to arrange the check this way, 
            j = negative_diag_stiff(i);                             % otherwise MATLAB gives problems of 
            vol.stiff(j,j) = abs(vol.stiff(j,j));                   % memory outage (8GB RAM)
        end
    end
    
    if not(dbg)
        clear negative_diag_stiff null_diag_stiff
    end
    
    if side==1
        rl = 'right';
    elseif side==2
        rl = 'left';
    else
        error("The considered side is neither 1 (right) or 2 (left). Check the calling to ea_genvat_aniso function is done properly")
    end
    save([options.root, options.patientname, filesep, 'stimulations', filesep, 'conduction_model_', rl, '.mat'], 'vol', 'cond', 'el_cond', '-v7.3')
end


%% ________________________________________________________________________
%% COMPUTE POTENTIAL AND ELECTRIC FIELD (ITS GRADIENT)
ea_dispt('Computing the potential based on stimulation...')

% ##########
% %% find boundary points in the volume of interest
% % TODO: include this section in the function to build the cylinder mesh (it should speed up a little)
% vol.boundary = unique(boundary(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3)));
% if dbg_vis
%     figure 
%     plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'b.')
%     hold on
%     plot3(vol.pos(vol.boundary, 1), vol.pos(vol.boundary, 2), vol.pos(vol.boundary, 3), 'r.')
%     title('boundary nodes')
%     legend('nodes inside the volume of interest','boundary nodes', 'location', 'northeast')
% end


%% define sources and contacts
if ~isfield(S, 'sources')
    S.sources=1:4;
end

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

% preallocate memory to store the gradient
gradient = cell(length(S.sources), 1);


%%
active_coords = [];
for source = S.sources  % ##### TODO: check if this block allows for multipolar stimulation, although it was not required by the project #####
    active_contacts = [];
    ix = [];
    voltix = [];
    stimsource = S.([sidec,'s',num2str(source)]);
    constvol = stimsource.va==1; % constvol is 1 for constant voltage and 0 for constant current.
    ac = find(S.activecontacts{side});
    U = zeros(length(ac),1);
    for con = ac  
        if S.([sidec,'s',num2str(source)]).amp && S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc  % find active contact corresponding to active source
            %% organize information to feed into the ea_apply_dbs function for potential computation
            
            % find the nodes of the active contact in scope
            if SI
                active_contacts = electrode.contacts(con).vertices / 1e3;  % in m
            else
                active_contacts = electrode.contacts(con).vertices;  % in mm
            end

            % find elements in mesh corresponding to nodes of the active
            % contact in scope
            vol_active = unique(knnsearch(vol.pos, active_contacts));
            
            % define the activeidx structure, that organizes the
            % information for stimulation in a way that fits ea_apply_dbs
            activeidx(source).con(con).ix = vol_active;
            activeidx(source).con(con).pol = S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).pol;
            activeidx(source).con(con).perc = S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc;
            
            if constvol
                U(con)=(logical(stimsource.(cnts{con}).perc))*stimsource.amp; % do not split amplitude in constant voltage setting.
            else
                U(con)=(stimsource.(cnts{con}).perc/100)*stimsource.amp;
            end

            if stimsource.(cnts{con}).pol==1
                U(con)=U(con)*-1;
            end
            
            if any(U>0) 
                unipolar=0;
                U=U/2;
            else
                unipolar=1;
            end
            
            active_coords = [active_coords; active_contacts(round(size(active_contacts,1)/3),:)]; % find coordinates where the contacts are active
%             active_coords = [active_coords; coords(logical(U), :) / 1e3];  % find coordinates where the contacts are active in m
%             active_coords = [active_coords; coords(logical(U), :)];  % find coordinates where the contacts are active in mm
            
            ix = [ix;activeidx(source).con(con).ix]; %#ok<*AGROW>  % belonging of an active element to a contact
            voltix_new = [repmat(U(con), length(activeidx(source).con(con).ix), 1), ...
                      repmat(con, length(activeidx(source).con(con).ix), 1)];   
            if ~constvol
                voltix_new(:,1) = voltix_new(:,1) / 1000;  % from mA to A
            end
            voltix = [voltix; voltix_new];  % value of stimulation of the active element
            
            if dbg_vis
                figure
                plot3(elmodel.ctr(:,1), elmodel.ctr(:,2), elmodel.ctr(:,3), 'b.')
                hold on
                plot3(vol.pos(ix,1), vol.pos(ix,2), vol.pos(ix,3), 'r.', 'MarkerSize', 20)
                plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'g.')
                legend('electrode', 'active', 'volume', 'location', 'northeast')
                title('active elements')
            end
        end
    end

    
    %% compute potential distribution and gradient
    if ~(isempty(active_contacts))          
        potential = ea_apply_dbs(vol,ix,voltix,unipolar,constvol,vol.boundary); % output in V. 4 indexes insulating material.
        gradient{source} = ea_calc_gradient(vol,potential);

% % ##### block below was commented #####
%         %% get high EF values for active electrodes
%         % this can be adjusted by assigning all the tetrahedra belonging to the
%         % active electrode a new value: gradient{source}(elec_tet_ix,:) = new_value;
%         elec_tet_ix = sub2ind(size(vol.pos), vertcat(ix,ix,ix), vertcat(ones(length(ix),1), ones(length(ix),1).*2, ones(length(ix),1).*3));
%         elec_tet_ix = find(sum(ismember(vol.tet,elec_tet_ix),2)==4);
% 
%         % gradient{source}(elec_tet_ix,:) = repmat(max(gradient{source}),[length(elec_tet_ix),1]); %assign maximum efield value
%         tmp = sort(abs(gradient{source}),'descend');
%         gradient{source}(elec_tet_ix,:) = repmat(mean(tmp(1:ceil(length(tmp(:,1))*0.001),:)),[length(elec_tet_ix),1]); % choose mean of highest 0.1% as new efield value
%         clear tmp

    else
        % if the source is not active, the solution is trivial
        gradient{source} = zeros(size(vol.tet,1),3);
    end
end

if not(dbg)
    clear voltix_new U activeidx elec_tet_ix ac con active_contacts voltix ix stimsource vol_active
end

% combine gradients from all sources
gradient=gradient{1}+gradient{2}+gradient{3}+gradient{4}; 

% ##########
% if SI
% % #####
% % convert back from m to mm
% vol.pos = vol.pos * 1e3;
% vol.ctr = vol.ctr * 1e3;
% vol.r = vol.r * 1e3;
% vol.unit = 'mm';
% % #####
% end


%% ________________________________________________________________________
%% FLOWFIELD VISUALIZATION
ea_dispt('Calculating quiver field of gradient for display purposes...');
midpts = vol.ctr;

if SI; midpts = midpts * 1e3; end  % convert back to mm

if dbg_vis
    figure
    quiver3(midpts(:,1),midpts(:,2),midpts(:,3),gradient(:,1),gradient(:,2),gradient(:,3))
    title('quiver field')
end


% % ##### block below added #####
%% remove electrode
if options.prefs.machine.vatsettings.aniso_removeElectrode
    gradient(el_cond,:) = 0;
    options.prefs.machine.vatsettings.aniso_removeElectrode = 0;
end


%% select subset of points to use for quiver representation
reduc=10;
indices = (1:reduc:length(midpts) + round(randn*(reduc/3)))';
indices=unique(indices(2:end-1));

% in order to have non-biased representation a component of randomness has
% been added. This requires the interventions of the following lines
indices(indices==0)=[];
indices(indices>length(midpts))=[];


%% transform to template space if necessary
if options.native==1 && options.orignative==0  % case if we are visualizing in MNI but want to calculate VTA in native space -> now transform back to MNI
    c=midpts';
    [~,anatpresent]=ea_assignpretra(options);
    V=ea_open_vol([options.root,options.patientname,filesep,anatpresent{1}]);  % this function obtains image volume information (header) of the files just found in the previous line
    c=V.mat\[c;ones(1,size(c,2))];  % V.mat is the pose of the image (wrt patient)
    midpts=ea_map_coords(c(1:3,:), ...  % Coordinates mapping
        [options.root,options.patientname,filesep,anatpresent{1}], ...
        [options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'], ...  % this file is generated by SPM
        '')';
    midpts=midpts(:,1:3);
    options.native=options.orignative;  % go back to template space
    clear c V
end

clear reduc

if dbg_vis && exist('c', 'var')
    figure
    plot3(midpts(:,1),midpts(:,2),midpts(:,3),'b.');
    title('Volume warped to native space')
end


%% ________________________________________________________________________
%% SAVE FUNCTION OUTPUTS
% no output to the user, since the following function is by itself very verbose
[vatfv,vatvolume,radius] = ea_write_vta_nii(S,stimname,midpts,indices,elspec,active_coords,voltix,constvol,thresh,vol,gradient,side,resultfig,options);


%% ________________________________________________________________________
%% OUTPUT
varargout{1}=vatfv;
varargout{2}=vatvolume;
varargout{3}=radius;
ea_dispt(''); % stop chain of timed processes.


%% 
if dbg
    fprintf("\nDEBUG ___________________________________________________\n")
    toc;
    fprintf("_________________________________________________________\n\n")
end


end  % function




function potential = ea_apply_dbs(vol,elec,val,unipolar,constvol,boundarynodes)
if constvol
    if unipolar
        dirinodes = [boundarynodes; elec];  % dirichlet nodes are the nodes for which the potential is fixed at dirival
    else
        dirinodes = elec;
    end

    rhs = zeros(length(vol.pos),1);  % right hand side vector
    dirival = zeros(size(vol.pos,1),1);
    dirival(elec) = val(:,1);
else

    if unipolar
        dirinodes = boundarynodes;
    else
        dirinodes = 1;
    end
    dirival = zeros(size(vol.pos,1),1);

    rhs = zeros(size(vol.pos,1),1);
    uvals=unique(val(:,2));
    if unipolar && length(uvals)==1
        elec_center_id = ea_find_elec_center(elec,vol.pos);
        rhs(elec_center_id) = val(1,1);
    else
        for v=1:length(uvals)
            elec_center_id = ea_find_elec_center(elec(val(:,2)==uvals(v)),vol.pos);
            thesevals=val(val(:,2)==uvals(v),1);
            rhs(elec_center_id) = thesevals(1);
        end

        %warning('Bipolar constant current stimulation currently not implemented!');
    end
end

% perform preconditioning on the stiffness matrix to make the next computation faster 
[stiff, rhs] = ea_dbs(vol.stiff,rhs,dirinodes,dirival);  

potential = ea_sb_solve(stiff,rhs);
end  % subfunction



% the following 3 functions are used in ea_apply_dbs in order to obtain the
% potential distribution
function center_id = ea_find_elec_center(elec, pos)

center = mean(pos(elec,:));

dist_center = sqrt(sum((pos(elec,:)-repmat(center,length(elec),1)).^2,2));
[~, elec_id] = min(dist_center);
center_id = elec(elec_id);
end  % subfunction



function [stiff,rhs] = ea_dbs(stiff,rhs,dirinodes,dirival)  % the function sets the boundary conditions making the siffness matrix informed of this

diagonal = diag(stiff);
stiff = stiff + stiff';
rhs = rhs - stiff*dirival;
stiff(dirinodes,:) = 0.0;
stiff(:,dirinodes) = 0.0;
diagonal = -diagonal;
%diagonal(:) = 0;
diagonal(dirinodes) = 1.0;
stiff = stiff + spdiags(diagonal(:),0,length(diagonal),length(diagonal));
rhs(dirinodes) = dirival(dirinodes);
end  % subfunction



function x = ea_sb_solve(sysmat,vecb) %#ok<INUSD>

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $
try
    L = ichol(sysmat); %#ok<NASGU>
catch
    alpha = max(sum(abs(sysmat),2)./diag(sysmat))-2;
    L = ichol(sysmat, struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %#ok<NASGU>
end

%scalen
[~,x]=evalc('pcg(sysmat,vecb,10e-10,5000,L,L'',vecb)');
end  % subfunction



% once I've computed the potential distribution with the 3 functions above, 
% I need to compute its gradient
function gradient = ea_calc_gradient(vol,potential)
normal = cross(vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:));
gradient = repmat(potential(vol.tet(:,1))./sum(normal.*(vol.pos(vol.tet(:,1),:)-(vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:),vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:));
gradient = gradient + repmat(potential(vol.tet(:,2))./sum(normal.*(vol.pos(vol.tet(:,2),:)-(vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:));
gradient = gradient + repmat(potential(vol.tet(:,3))./sum(normal.*(vol.pos(vol.tet(:,3),:)-(vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:));
gradient = gradient + repmat(potential(vol.tet(:,4))./sum(normal.*(vol.pos(vol.tet(:,4),:)-(vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:))/3),2),1,3).*normal;
end  % subfunction

