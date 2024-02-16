      
function chanMapFilePath = JF_imroToChannelMapLoc(channelMapIMRO, metaFile, recompute )
if recompute
    metaInfo = dir(metaFile);
    chanMapFilePath = SGLXMetaToCoords_JF(metaInfo.name, metaInfo.folder, metaInfo.folder);
else
    if contains(channelMapIMRO, 'NPtype24_hStripe_')
        chanMap = "chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype21')
        chanMap = "chanMapNP2_1Shank_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_hStripe_shanks01')
        chanMap = "chanMapNP2_4Shank_bottRow_shank01_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_hStripe_shanks23') 
        chanMap = "chanMapNP2_4Shank_bottRow_shank23_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_shank0')
        chanMap = "chanMapNP2_4Shank_bottRow_shank0_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_shank01') 
        chanMap = "chanMapNP2_4Shank_bottRow_shank1_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_shank2')
        chanMap = "chanMapNP2_4Shank_bottRow_shank2_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype24_shank3') 
        chanMap = "chanMapNP2_4Shank_bottRow_shank3_flipper0x2Emat_kilosortChanMap.mat";
    elseif contains(channelMapIMRO, 'NPtype3B_')
        chanMap = "chanMapNP1_bottRow_flipper.mat" ;
    else
        error('channel map not recognized')
    end
    chanMapFilePath = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles', chanMapFile];

end
end
