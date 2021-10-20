%load probe_ccf for each animal, find which contains regions and at what
%depth 

mice = { 'JF054'}; % JF020 doesn't have histology yet. JF021 and JF022 to come, and many more
locations = {'CP', 'GPe', 'SNr'};

allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF';
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
clearvars regionID
for iRegion = 1:size(locations,2)
    regionID(iRegion) = find(strcmp(st.acronym, locations{iRegion}));
end
for iMouse = [1:16]
    %load probe_ccf: regions and depths 
    clearvars depthBorders
    animal = mice{iMouse};
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)
    for iProbe = [1:size(probe_ccf,1)]
        for iRegion = 1:size(locations,2)
            isPresent = find(probe_ccf(iProbe).trajectory_areas == regionID(iRegion));
            if ~isempty(isPresent)
                depthBorders(:,iProbe,iRegion) = [probe_ccf(iProbe).probe_depths(isPresent(1)), ...
                    probe_ccf(iProbe).probe_depths(isPresent(end))];
            end
        end
    end
end
