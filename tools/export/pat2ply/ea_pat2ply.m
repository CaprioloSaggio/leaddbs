function ea_pat2ply(uipatdir,handles)




    
atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};


cfv(1)=ea_atlas2ply({atlasset},[uipatdir,filesep,'export',filesep,'ply',filesep,'anatomy.ply']);
try % this is DBS specific.
    options=ea_detsides(ea_getptopts(uipatdir));
    cnt=1;
    for side=options.sides
        cfv(1+cnt)=ea_electrode2ply([uipatdir,filesep],side,handles);
    end
    cfvel=ea_concatfv(cfv(2:end));
    plywrite([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_electrodes.ply'],cfvel.faces,cfvel.vertices,cfvel.facevertexcdata,repmat(100,size(cfvel.facevertexcdata,1),1));
    
    cfv=ea_concatfv(cfv);
    plywrite([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.faces,cfv.vertices,cfv.facevertexcdata,repmat(100,size(cfv.facevertexcdata,1),1));
    %write_ply(cfv.vertices',cfv.faces',[uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply']);
    %ea_patch2ply([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.vertices',cfv.faces',cfv.facevertexcdata');
end


