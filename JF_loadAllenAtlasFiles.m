function [tv, av, st, bregma] = JF_loadAllenAtlasFiles(probe2ephysFile)
myPaths;
if contains(probe2ephysFile, 'brainreg') % new brainreg registration

    atlasResolution_um = 25; % voxel size. currently in all dimensions, both for atlas and acquired image
    atlasSpecies = 'mouse'; % atlas species
    atlasType = 'allen'; % atlas name

    atlasLocation = dir([brainglobeLocation, atlasType, '_', ...
        atlasSpecies, '_', num2str(atlasResolution_um), 'um*']); % atlas location
    [tv, av, st, bregma] = bd_loadAllenAtlas([atlasLocation.folder, filesep, atlasLocation.name]);
else % old elastix-based registration

    allen_atlas_path = [allenAtlasPath, filesep, 'allenCCF/'];
    tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
    av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
    error('JF_loadAllData not setup for elastix data yet')
end
end