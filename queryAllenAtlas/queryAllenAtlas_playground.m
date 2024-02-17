%docs: https://allensdk.readthedocs.io/en/latest/unionizes.html

%% info
cl_myPaths;
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;

curr_plot_structure = st.id(find(contains(st.name, 'audoputamen')));

%% dms and ps exp ids
expIDs_medialStriatum = [160537796, ... %anterior
    124059700,...
    159223001, ... %anterior
    146553266, ... %anterior
    112458831,...
    100142580,...
    301620241,...
    278434443,...
    127140981,...
    308395312,...
    278434443,...
    159941339,...
    307910595,...
    112307754,...
    293366741,...
    293366035,...
    160540013,...
    301180385,...
    112458831, ... 
    127762867, ... %posterior
    117317884, ... %posterior
    120762196]; %posterior

%% get all VIS Ctx experiments
injectionAreas = {'VISp', 'VISl', 'VISpl', 'VISpm', 'VISal', 'VISam', 'CP'};
injectionColors = [rgb('DarkRed'); rgb('OrangeRed'); rgb('DarkGreen'); rgb('MediumBlue'); rgb('Purple'); rgb('HotPink'); rgb('SkyBlue')];
viewVals = [0, 90; 0, 0; 90, 0];
viewLegends = {'AP (\mum)', 'ML (\mum)'; 'AP (\mum)', 'DV (\mum)'; 'ML (\mum)', 'DV (\mum)'};
collapseHemispheres = 1; %plot only one hemispehre
cmapData = brewermap(100, 'Reds');

%% save data
visProjectionData = struct;
visProjectionData.max_voxel_x = [];
visProjectionData.max_voxel_y = [];
visProjectionData.max_voxel_z = [];
visProjectionData.normalized_projection_volume = [];
visProjectionData.projection_intensity = [];
visProjectionData.injection_area = [];

cp_gpe_ProjectionData = struct;
cp_gpe_ProjectionData.max_voxel_x = [];
cp_gpe_ProjectionData.max_voxel_y = [];
cp_gpe_ProjectionData.max_voxel_z = [];
cp_gpe_ProjectionData.normalized_projection_volume = [];
cp_gpe_ProjectionData.projection_intensity = [];
cp_gpe_ProjectionData.injection_area = [];

cp_snr_ProjectionData = struct;
cp_snr_ProjectionData.max_voxel_x = [];
cp_snr_ProjectionData.max_voxel_y = [];
cp_snr_ProjectionData.max_voxel_z = [];
cp_snr_ProjectionData.normalized_projection_volume = [];
cp_snr_ProjectionData.projection_intensity = [];
cp_snr_ProjectionData.injection_area = [];

for iInjection = 1:size(injectionAreas, 2)


    if ~exist([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'], 'file')
        if strcmp(injectionAreas(iInjection), 'CP')
            expIDs = expIDs_medialStriatum;
        else
            expIDs = findAllenExperiments('injection', injectionAreas{iInjection}, 'primary', true);
        end

        proj = [];
        for iExpID = 1:size(expIDs, 2)
            proj_temp = getProjectionDataFromExperiment(expIDs(iExpID));
            proj = [proj, proj_temp{1, 1}];
        end
        save([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/DATA/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'], 'proj')
    else
        load([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/DATA/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'])
    end

    if strcmp(injectionAreas(iInjection), 'CP')
        curr_plot_structure = st.id(find(contains(st.acronym, 'GPe')));
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;
    
        cp_gpe_ProjectionData.max_voxel_x = [cp_gpe_ProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        cp_gpe_ProjectionData.max_voxel_y = [cp_gpe_ProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        cp_gpe_ProjectionData.max_voxel_z = [cp_gpe_ProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        cp_gpe_ProjectionData.normalized_projection_volume = [cp_gpe_ProjectionData.normalized_projection_volume, [proj(theseTargets).normalized_projection_volume]];
        cp_gpe_ProjectionData.projection_intensity = [cp_gpe_ProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        cp_gpe_ProjectionData.injection_area = [cp_gpe_ProjectionData.injection_area, areaName];
        
        curr_plot_structure = st.id(find(contains(st.acronym, 'SNr')));
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;
    
        cp_snr_ProjectionData.max_voxel_x = [cp_snr_ProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        cp_snr_ProjectionData.max_voxel_y = [cp_snr_ProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        cp_snr_ProjectionData.max_voxel_z = [cp_snr_ProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        cp_snr_ProjectionData.normalized_projection_volume = [cp_snr_ProjectionData.normalized_projection_volume, [proj(theseTargets).normalized_projection_volume]];
        cp_snr_ProjectionData.projection_intensity = [cp_snr_ProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        cp_snr_ProjectionData.injection_area = [cp_snr_ProjectionData.injection_area, areaName];
    else
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;
    
        visProjectionData.max_voxel_x = [visProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        visProjectionData.max_voxel_y = [visProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        visProjectionData.max_voxel_z = [visProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        visProjectionData.normalized_projection_volume = [visProjectionData.normalized_projection_volume, [proj(theseTargets).normalized_projection_volume]];
        visProjectionData.projection_intensity = [visProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        visProjectionData.injection_area = [visProjectionData.injection_area, areaName];
    end

    clearvars areaName proj


end
save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/VIS_Str.mat'], 'visProjectionData')
save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/Str_gpe_medial.mat'], 'cp_gpe_ProjectionData')
save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/Str_snr_medial.mat'], 'cp_snr_ProjectionData')
%% plot data 3D - striatum 
curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;

curr_plot_structure = st.id(find(contains(st.name, 'audoputamen')));

cF=figure();

for iInjection = 1:size(injectionAreas, 2)-1

    load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/VIS_Str.mat'])

    for iView = 1%:3
        h1 = subplot(1, 1, iView);
        if iInjection==1
        [~, brain_outline] = plotBrainGrid([], h1); hold on;
            structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                'Faces', structure_3d.faces, ...
                'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
     
        % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
        %            'Faces', structure_3d.faces, ...
        %            'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
        end
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])

        hold on;
        axis equal
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        %lim([-10, ap_max + 10])
        % if collapseHemispheres
        %     ylim([-10, ml_max / 2 + 10])
        % else
        %     ylim([-10, ml_max + 10])
        % end

        zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        %         if iInjection == 1
        %             xlabel(viewLegends{iView, 1})
        %             ylabel(viewLegends{iView, 2})
        %         end

        % xticks([0 5 10])
        % xticklabels({'x = 0','x = 5','x = 10'})
        projData = visProjectionData;
        if collapseHemispheres
            max_voxel_x = [projData.max_voxel_x] / 10;

            max_voxel_z = [projData.max_voxel_z] / 10;
            z_rel = max_voxel_z - ml_max / 2;
            max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);

            max_voxel_y = [projData.max_voxel_y] / 10;
        else
            max_voxel_x = [projData.max_voxel_x] / 10;
            max_voxel_y = [projData.max_voxel_y] / 10;
            max_voxel_z = [projData.max_voxel_z] / 10;
        end
        %[~, brain_outline] = plotBrainGrid([],[]);
                % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                %     'Faces', structure_3d.faces, ...
                %     'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])
                % 

        theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        transp = [projData.projection_intensity(theseTargets)] ./ (max([projData.projection_intensity(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt == 0) = 1;
        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [projData.normalized_projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
            1)
        if iView ~= 1
            set(h1, 'Ydir', 'reverse')
        end
        set(h1, 'box', 'off', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', [])
        %xlim([300, 800])
        set(gca, 'Visible', 'off')
        makepretty;
        %title(injectionAreas{iInjection})
    end
end

view([0, 0])
%save as .avi rotating vid
set(gcf, 'Color', 'w')
OptionZ.FrameRate = 25;
OptionZ.Duration = 15;
OptionZ.Periodic = true;
CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_vis_to_str'], OptionZ)

%% plot data 3D - GPe
figure();
curr_plot_structure_idx = find(contains(st.acronym, 'GPe'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;

curr_plot_structure = st.id(find(contains(st.acronym, 'GPe')));

for iInjection = size(injectionAreas, 2)

    load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/Str_gpe_medial.mat'])

    for iView = 1%:3
        h2 = subplot(1, 1, iView);
        if iInjection==size(injectionAreas, 2)
        [~, brain_outline] = plotBrainGrid([], h2); hold on;
              structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                'Faces', structure_3d.faces, ...
                'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
     
        % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
        %            'Faces', structure_3d.faces, ...
        %            'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
        end
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])

        hold on;
        axis equal
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        %lim([-10, ap_max + 10])
        % if collapseHemispheres
        %     ylim([-10, ml_max / 2 + 10])
        % else
        %     ylim([-10, ml_max + 10])
        % end

        zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        %         if iInjection == 1
        %             xlabel(viewLegends{iView, 1})
        %             ylabel(viewLegends{iView, 2})
        %         end

        % xticks([0 5 10])
        % xticklabels({'x = 0','x = 5','x = 10'})
        projData = cp_gpe_ProjectionData;
        if collapseHemispheres
            max_voxel_x = [projData.max_voxel_x] / 10;

            max_voxel_z = [projData.max_voxel_z] / 10;
            z_rel = max_voxel_z - ml_max / 2;
            max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);

            max_voxel_y = [projData.max_voxel_y] / 10;
        else
            max_voxel_x = [projData.max_voxel_x] / 10;
            max_voxel_y = [projData.max_voxel_y] / 10;
            max_voxel_z = [projData.max_voxel_z] / 10;
        end
        %[~, brain_outline] = plotBrainGrid([],[]);
                % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                %     'Faces', structure_3d.faces, ...
                %     'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])
                % 

        theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        transp = [projData.projection_intensity(theseTargets)] ./ (max([projData.projection_intensity(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt == 0) = 1;
        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [projData.normalized_projection_volume(theseTargets)]*100, cmapData(round(round(transp, 2).*100), :), 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
            1)
        if iView ~= 1
            set(h2, 'Ydir', 'reverse')
        end
        set(h2, 'box', 'off', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', [])
        %xlim([300, 800])
        set(gca, 'Visible', 'off')
        makepretty;
        %title(injectionAreas{iInjection})
    end
end

view([0, 0])
%save as .avi rotating vid
set(gcf, 'Color', 'w')
OptionZ.FrameRate = 25;
OptionZ.Duration = 15;
OptionZ.Periodic = true;
CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_str_to_gpe'], OptionZ)

%% plot data 3D -SNr

figure();
curr_plot_structure_idx = find(contains(st.acronym, 'SNr'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;

curr_plot_structure = st.id(find(contains(st.acronym, 'SNr')));

for iInjection = size(injectionAreas, 2)

    load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/Str_snr_medial.mat'])

    for iView = 1%:3
        h2 = subplot(1, 1, iView);
        if iInjection==size(injectionAreas, 2)
        [~, brain_outline] = plotBrainGrid([], h2); hold on;
              structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                'Faces', structure_3d.faces, ...
                'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
     
        % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
        %            'Faces', structure_3d.faces, ...
        %            'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
        end
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])

        hold on;
        axis equal
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        %lim([-10, ap_max + 10])
        % if collapseHemispheres
        %     ylim([-10, ml_max / 2 + 10])
        % else
        %     ylim([-10, ml_max + 10])
        % end

        zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        %         if iInjection == 1
        %             xlabel(viewLegends{iView, 1})
        %             ylabel(viewLegends{iView, 2})
        %         end

        % xticks([0 5 10])
        % xticklabels({'x = 0','x = 5','x = 10'})
        projData = cp_snr_ProjectionData;
        if collapseHemispheres
            max_voxel_x = [projData.max_voxel_x] / 10;

            max_voxel_z = [projData.max_voxel_z] / 10;
            z_rel = max_voxel_z - ml_max / 2;
            max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);

            max_voxel_y = [projData.max_voxel_y] / 10;
        else
            max_voxel_x = [projData.max_voxel_x] / 10;
            max_voxel_y = [projData.max_voxel_y] / 10;
            max_voxel_z = [projData.max_voxel_z] / 10;
        end
        %[~, brain_outline] = plotBrainGrid([],[]);
                % structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
                %     'Faces', structure_3d.faces, ...
                %     'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
                % minV = min(min(structure_patch.Vertices(:,:)));
                % maxV = max(max(structure_patch.Vertices(:,:)));
                % xlim([])
                % 

        theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        transp = [projData.projection_intensity(theseTargets)] ./ (max([projData.projection_intensity(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt == 0) = 1;
        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [projData.normalized_projection_volume(theseTargets)]*100, cmapData(round(round(transp, 2).*100), :), 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
            1)
        if iView ~= 1
            set(h2, 'Ydir', 'reverse')
        end
        set(h2, 'box', 'off', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', [])
        %xlim([300, 800])
        set(gca, 'Visible', 'off')
        makepretty;
        %title(injectionAreas{iInjection})
    end
end


view([0, 0])
%save as .avi rotating vid
set(gcf, 'Color', 'w')
OptionZ.FrameRate = 25;
OptionZ.Duration = 15;
OptionZ.Periodic = true;
CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_str_to_snr'], OptionZ)


figure();
subplot(144)
colormap(cmapData)
h = colorbar;
h.Title.String = ['Normalized', newline, 'projection', newline, 'intensity'];
makepretty;
set(gcf, 'Color', 'white')

