allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure = st.id(find(contains(st.name, 'audoputamen')));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
structure_alpha = 0.2;

%get all VIS Ctx experiments
injectionAreas = {'VISp', 'VISl',  'VISpl', 'VISpm', 'VISal', 'VISam'};
injectionColors = [rgb('DarkRed'); rgb('OrangeRed'); rgb('DarkGreen'); rgb('MediumBlue'); rgb('Purple'); rgb('HotPink'); rgb('SkyBlue')];
viewVals = [0,90;0,0;90,0];
figure();
for iInjection = 1:size(injectionAreas,2)
    expIDs = findAllenExperiments('injection', injectionAreas{iInjection}, 'line', '0', 'primary', true);
    
    proj = getProjectionDataFromExperiment(expIDs(1));

    for iView = 1:3
        ss = subplot(3, size(injectionAreas, 2), (size(injectionAreas, 2))*(iView-1)+iInjection);
        %[~, brain_outline] = plotBrainGrid([], ss);
        hold on;
        view(viewVals(iView,:))
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        ylim([-10, ml_max + 10])
        zlim([-10, dv_max + 10])
        structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
            'Faces', structure_3d.faces, ...
            'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);


        projData = proj{1, 1};
        theseTargets = ismember([projData.structure_id], double(curr_plot_structure)) & [projData.max_voxel_density] > 0;
        hold on;
        scatter3([projData(theseTargets).max_voxel_x]/10, [projData(theseTargets).max_voxel_z]/10, [projData(theseTargets).max_voxel_y]/10, ...
            [projData(theseTargets).normalized_projection_volume]*10, injectionColors(iInjection,:), 'filled')
    end
end