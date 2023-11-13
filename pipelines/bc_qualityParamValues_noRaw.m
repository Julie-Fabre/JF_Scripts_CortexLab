function param = bc_qualityParamValues_noRaw(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV)
    param = bc_qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV);
    param.removeDuplicateSpikes = 0; 
    param.extractRaw = 0;
end
