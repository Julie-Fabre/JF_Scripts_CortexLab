animalsType = {'Naive'};
regionsNames = {'CP', 'CP', 'STN', 'GPe', 'SNr', 'GPi'};
regions = {'DMS', 'PS', 'STN', 'GPe', 'SNr', 'GPi'};
recordingInfo = readtable('C:\Users\Julie\Dropbox\Analysis\Recordings - Sheet1.csv');


for iType = 1:size(animalsType, 2)
    theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    theseColors = {rgb('DeepSkyBlue'); rgb('Turquoise'); rgb('MediumOrchid'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink')};

    allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
    tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
    av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
    slice_spacing = 10;
    structure_alpha = 0.2;
    %get colors (overrride allen)

    figure();
    [~, brain_outline] = plotBrainGrid([], []);

    %overlay regions
    for iRegion = 1:size(regions, 2)

        curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);

        hold on;
        axis vis3d equal off manual
        view([-30, 25]);
        caxis([0, 300]);
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        ylim([-10, ml_max + 10])
        zlim([-10, dv_max + 10])
        structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
            'Faces', structure_3d.faces, ...
            'FaceColor', theseColors{iRegion, :}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
        %plot probe tracks
        theseProbes = strcmp(recordingInfo.Location, regions{iRegion});
        %get animal and probe, load track
        
    end
end
