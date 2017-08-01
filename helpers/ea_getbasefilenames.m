function [basefilenames,mostallfiles]=ea_getbasefilenames(options,append)

if exist('options','var')
    if ~isempty(options)
        prefs=ea_prefs(options.patientname);
    else
        prefs=ea_prefs('');
    end
else
    prefs=ea_prefs('');
end


basefilenames={
    prefs.prenii_unnormalized,...
    prefs.prenii_unnormalized_t1,...
    prefs.prenii_unnormalized_pd,...
    prefs.tranii_unnormalized,...
    prefs.cornii_unnormalized,...
    prefs.sagnii_unnormalized,...
    prefs.rawctnii_unnormalized,...
    prefs.ctnii_coregistered,...
    prefs.rest_default,...
    prefs.dti
};

mostallfiles={
    prefs.prenii,...
    prefs.tranii,...
    prefs.cornii,...
    prefs.sagnii,...
    prefs.ctnii,...
    prefs.gprenii,...
    prefs.gtranii,...
    prefs.gcornii,...
    prefs.gsagnii,...
    prefs.gctnii,...
    prefs.tp_ctnii_coregistered,...
    prefs.tp_gctnii,...
    prefs.b0,...
    prefs.fa,...
    prefs.fa2anat,...
    ['l', prefs.fa2anat],...
    ['gl', prefs.fa2anat],...
    'grid.nii',...
    'glgrid.nii',...
};

if exist('append','var')
   if append
    for e=1:length(basefilenames)
       thisfi=basefilenames{e};
       [~,base,ext]=fileparts(thisfi);
       basefilenames{e}=[base,num2str(append),ext];
    end
    for e=1:length(mostallfiles)
        thisfi=mostallfiles{e};
        [~,base,ext]=fileparts(thisfi);
        mostallfiles{e}=[base,num2str(append),ext];
    end
   end
end
