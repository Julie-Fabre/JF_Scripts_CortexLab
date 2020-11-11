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
injectionAreas = {'VISp', 'VISl', 'VISpl', 'VISpm', 'VISal', 'VISam'};
injectionColors = [rgb('DarkRed'); rgb('OrangeRed'); rgb('DarkGreen'); rgb('MediumBlue'); rgb('Purple'); rgb('HotPink'); rgb('SkyBlue')];
viewVals = [0, 90; 0, 0; 90, 0];
viewLegends = {'AP (\mum)', 'ML (\mum)'; 'AP (\mum)', 'DV (\mum)'; 'ML (\mum)', 'DV (\mum)'};
collapseHemispheres = 0; %plot only one hemispehre
figure();
for iInjection = 1:size(injectionAreas, 2)
    if ~exist(['C:\Users\Julie\Dropbox\MATLAB\JF_scripts_CortexLab\queryAllenAtlas\', injectionAreas{iInjection}, '.mat'], 'file')
        expIDs = findAllenExperiments('injection', injectionAreas{iInjection}, 'primary', true);

        proj = [];
        for iExpID = 1:size(expIDs, 2)
            proj_temp = getProjectionDataFromExperiment(expIDs(iExpID));
            proj = [proj, proj_temp{1, 1}];
        end
        save(['C:\Users\Julie\Dropbox\MATLAB\JF_scripts_CortexLab\queryAllenAtlas\', injectionAreas{iInjection}, '.mat'], 'proj')
    else
        load(['C:\Users\Julie\Dropbox\MATLAB\JF_scripts_CortexLab\queryAllenAtlas\', injectionAreas{iInjection}, '.mat'])
    end
    for iView = 1:3
        ss = subplot(3, size(injectionAreas, 2), (size(injectionAreas, 2))*(iView - 1)+iInjection);
        %[~, brain_outline] = plotBrainGrid([], ss);
        hold on;
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        ylim([-10, ml_max + 10])
        zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        if iInjection == 1
            xlabel(viewLegends{iView, 1})
            ylabel(viewLegends{iView, 2})
        end
        projData = proj;
        if collapseHemispheres
            max_voxel_x = [projData.max_voxel_x]/10;
           
            max_voxel_z = [projData.max_voxel_z]/10;
            z_rel = max_voxel_z - ml_max / 2;
            max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);
            
            max_voxel_y = [projData.max_voxel_y]/10;
        else
            max_voxel_x = [projData.max_voxel_x]/10;
            max_voxel_y = [projData.max_voxel_y]/10;
            max_voxel_z = [projData.max_voxel_z]/10;
        end

        structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
            'Faces', structure_3d.faces, ...
            'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);


        
        theseTargets = ismember([projData.structure_id], double(curr_plot_structure)) & [projData.max_voxel_density] > 0;
        hold on;

        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [projData(theseTargets).normalized_projection_volume]*10, injectionColors(iInjection, :), 'filled')
    end
end