function [success,commands]=ea_checkinstall(cmd,force,checkonly,robot)
success=1;
if ~exist('checkonly','var')
    checkonly=0;
end
if ~exist('robot','var')
    robot=0;
end
earoot=ea_getearoot;
if ~exist('force','var')
    force=0;
end
switch cmd
    case 'list' % simply return list of installable datasets
        success={'Big Brain 100um subcortical (Amunts 2013)','Lead-DBS Macaque Toolbox','Structural group connectome (Horn 2013)'};
        commands={'bigbrain','macaque','groupconnectome2013'};
    case 'bigbrain'
        checkf=[earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly;
            success=~force;
            return;
        end
        if force==-1;
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('Bigbrain 100um subcortical',...
                [earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz'],...
                'bigbrain');
        else
            disp('BigBrain is installed.')
        end
    case 'macaque'
        checkf=[earoot,'toolbox',filesep,'macaque'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly;
            success=~force;
            return;
        end
        if force==-1;
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('Lead-DBS Macaque toolbox',...
                [earoot,'toolbox',filesep,'macaque.zip'],...
                'macaque');
        else
            disp('Macaque toolbox is installed.')
        end
    case 'groupconnectome2013'
        checkf=[ea_getconnectomebase('dmri'),'Groupconnectome (Horn 2013) full.mat'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly;
            success=~force;
            return;
        end
        if force==-1;
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('structural group connectome (Horn 2013)',...
                [ea_getconnectomebase('dmri'),'groupconnectome2013.zip'],...
                'group');
        else
            disp('Group Connectome (Horn 2013) is installed.')
        end
    otherwise
        success=0;
end


function success=ea_downloadasset(assetname,destination,id)
downloadurl = 'http://www.lead-dbs.org/release/download.php';
success=1;
disp(['Downloading ',assetname,'...'])
try
    webopts=weboptions('Timeout',5);
    websave(destination,downloadurl,'id',id,webopts);
catch
    try
        urlwrite([downloadurl,'?id=',id],destination,'Timeout',5);
    catch
        success=0;
    end
end

if success
    [~,~,ext] = fileparts(destination);
    if strcmp(ext,'.gz')
        gunzip(destination);
    elseif strcmp(ext,'.zip')
        unzip(destination);
    end
end

ea_delete(destination);


function force=ea_alreadyinstalled(checkf,checkonly,robot)
if ~exist(checkf,'file') % then file not there, should install anyways.
    force=1;
    return
end
if checkonly % return here.
    force=0;
    return
end
if ~robot
choice = questdlg('This dataset seems already to be installed. Do you wish to re-download it?', ...
    'Redownload Dataset?', ...
    'Yes','No','No');
if strcmp(choice,'Yes')
    force=1;
else
    force=-1;
end
else
    force=0;
end
