%docs: https://allensdk.readthedocs.io/en/latest/unionizes.html
function cl_getAllenAtlasProjections_visual(regionType, plotType, plotInjections)
% regionType : VISp, VIS, LGd, SCs,
% all_naive, mPFC, all
if nargin < 3 || isempty(plotInjections)
    plotInjections = 0;
end
if nargin < 2 || isempty(plotType)
    plotType = '2D';
end
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
expIDs_medialStriatum = [ ...
    124059700, ...
    159223001, ... %anterior
    146553266, ... %keep? anterior
    112458831, ... %keep?
    100142580, ... %keep??
    301620241, ... %keep?
    278434443, ...
    127140981, ... %keep??
    308395312, ... %keep??
    159941339, ...
    307910595, ... %keep??
    112307754, ... %keep?
    293366741, ...
    293366035, ...
    160540013, ...
    301180385, ...
    112458831, ...
    127762867, ... %posterior
    117317884, ... %posterior
    120762196]; %posterior

expIDs_centralStriatum = [307910595, ...
    272697238, ...
    264095536, ...
    112307754, ...
    120570964, ...
    160537796, ...
    127140981, ...
    157911832];


expIDs_medialPFC = [183618139, ...
    112952510, ...
    177459319, ...
    184169615, ...
    166460484, ...
    166055636, ...
    182814071, ...
    292172846, ...
    298000163];

%% get all VIS Ctx experiments
if strcmp(regionType, 'VISp')
    injectionAreas = {'VISp', 'DMS'};
elseif strcmp(regionType, 'VIS')
    injectionAreas = {'VISp', 'VISl', 'VISal', 'VISpm', 'VISam', 'DMS'}; %exclud VISp and others  for now : known small + axon go through posterior
elseif strcmp(regionType, 'LGd')
    injectionAreas = {'LGd', 'DMS'};
elseif strcmp(regionType, 'SCs')
    injectionAreas = {'SCs', 'DMS'};
elseif strcmp(regionType, 'mPFC')
    injectionAreas = {'MOs', 'DCS'};
elseif strcmp(regionType, 'all_naive')
    injectionAreas = {'VISpm', 'VISam', 'LGd', 'SCs', 'DMS'};
elseif strcmp(regionType, 'all')
    injectionAreas = {'VISpm', 'VISam', 'LGd', 'SCs', 'MOs', 'DS'};
else
    if iscell(regionType)
        injectionAreas = regionType;
    else
        injectionAreas = {regionType, 'DMS'};
    end
end

viewVals = [0, 90; 0, 0; 90, 0];
collapseHemispheres = 1; %plot only one hemispehre
cmapData = brewermap(100, 'Reds');

%% save data

visProjectionData = struct;
visProjectionData.max_voxel_x = [];
visProjectionData.max_voxel_y = [];
visProjectionData.max_voxel_z = [];
visProjectionData.projection_volume = [];
visProjectionData.volume = [];
visProjectionData.projection_intensity = [];
visProjectionData.injection_area = [];

cp_gpe_ProjectionData = struct;
cp_gpe_ProjectionData.max_voxel_x = [];
cp_gpe_ProjectionData.max_voxel_y = [];
cp_gpe_ProjectionData.max_voxel_z = [];
cp_gpe_ProjectionData.projection_volume = [];
cp_gpe_ProjectionData.volume = [];
cp_gpe_ProjectionData.projection_intensity = [];
cp_gpe_ProjectionData.injection_area = [];

cp_snr_ProjectionData = struct;
cp_snr_ProjectionData.max_voxel_x = [];
cp_snr_ProjectionData.max_voxel_y = [];
cp_snr_ProjectionData.max_voxel_z = [];
cp_snr_ProjectionData.projection_volume = [];
cp_snr_ProjectionData.volume = [];
cp_snr_ProjectionData.projection_intensity = [];
cp_snr_ProjectionData.injection_area = [];

for iInjection = 1:size(injectionAreas, 2)


    if ~exist([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/', injectionAreas{iInjection}, '.mat'], 'file')
        if strcmp(injectionAreas(iInjection), 'DMS')
            expIDs = expIDs_medialStriatum;
        elseif strcmp(injectionAreas(iInjection), 'DMS_single') %sanity check
            expIDs = expIDs_medialStriatum(1);
        elseif strcmp(injectionAreas(iInjection), 'DCS')
            expIDs = expIDs_centralStriatum;
        elseif strcmp(injectionAreas(iInjection), 'DS')
            expIDs = [expIDs_medialStriatum, expIDs_centralStriatum];
        elseif strcmp(injectionAreas(iInjection), 'MOs')
            expIDs = expIDs_medialPFC;
        else
            expIDs = findAllenExperiments('injection', injectionAreas{iInjection}, 'primary', true);
        end

        proj = [];
        for iExpID = 1:size(expIDs, 2)
            proj_temp = getProjectionDataFromExperiment(expIDs(iExpID));
            proj = [proj, proj_temp{1, 1}];
        end
        save([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/', injectionAreas{iInjection}, '.mat'], 'proj')
    else
        load([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/', injectionAreas{iInjection}, '.mat'])
    end
    % projection density: number of projecting pixels / voxel volume
    % injection density: number of projecting pixels in injection site / voxel volume
    if ismember(injectionAreas(iInjection), {'DMS', 'DCS', 'DS', 'DMS_single'})
        curr_plot_structure = st.id(find(contains(st.acronym, 'GPe')));
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;

        cp_gpe_ProjectionData.max_voxel_x = [cp_gpe_ProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        cp_gpe_ProjectionData.max_voxel_y = [cp_gpe_ProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        cp_gpe_ProjectionData.max_voxel_z = [cp_gpe_ProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        cp_gpe_ProjectionData.volume = [cp_gpe_ProjectionData.volume, [proj(theseTargets).volume]];

        cp_gpe_ProjectionData.projection_volume = [cp_gpe_ProjectionData.projection_volume, [proj(theseTargets).projection_volume]];
        cp_gpe_ProjectionData.projection_intensity = [cp_gpe_ProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        cp_gpe_ProjectionData.injection_area = [cp_gpe_ProjectionData.injection_area, areaName];

        curr_plot_structure = st.id(find(contains(st.acronym, 'SNr')));
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;

        cp_snr_ProjectionData.max_voxel_x = [cp_snr_ProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        cp_snr_ProjectionData.max_voxel_y = [cp_snr_ProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        cp_snr_ProjectionData.max_voxel_z = [cp_snr_ProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        cp_snr_ProjectionData.projection_volume = [cp_snr_ProjectionData.projection_volume, [proj(theseTargets).projection_volume]];
        cp_snr_ProjectionData.volume = [cp_snr_ProjectionData.volume, [proj(theseTargets).volume]];
        cp_snr_ProjectionData.projection_intensity = [cp_snr_ProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        cp_snr_ProjectionData.injection_area = [cp_snr_ProjectionData.injection_area, areaName];
    else
        theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;

        visProjectionData.max_voxel_x = [visProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
        visProjectionData.max_voxel_y = [visProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
        visProjectionData.max_voxel_z = [visProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
        visProjectionData.projection_volume = [visProjectionData.projection_volume, [proj(theseTargets).projection_volume]];
        visProjectionData.volume = [visProjectionData.volume, [proj(theseTargets).volume]];
        visProjectionData.projection_intensity = [visProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
        areaName = cell(1, size(find(theseTargets), 2));
        areaName(:) = {injectionAreas{iInjection}};
        visProjectionData.injection_area = [visProjectionData.injection_area, areaName];
    end

    clearvars areaName proj


end

%save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/VIS_Str.mat'], 'visProjectionData')
%save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/Str_gpe_medial.mat'], 'cp_gpe_ProjectionData')
%save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/Str_snr_medial.mat'], 'cp_snr_ProjectionData')
if strcmp(plotType, '3D')

    %% plot data 3D - striatum
    curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
    slice_spacing = 10;
    plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

    structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
        1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
    structure_alpha = 0.2;

    curr_plot_structure = st.id(find(contains(st.name, 'audoputamen')));

    cF = figure();

    for iInjection = 1:size(injectionAreas, 2) - 1

        %load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/VIS_Str.mat'])

        for iView = 1 %:3
            h1 = subplot(1, 1, iView);
            if iInjection == 1
                [~, brain_outline] = plotBrainGrid([], h1);
                hold on;
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
                [projData.projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
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
    %CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_vis_to_str'], OptionZ)

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

        %load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/Str_gpe_medial.mat'])

        for iView = 1 %:3
            h2 = subplot(1, 1, iView);
            if iInjection == size(injectionAreas, 2)
                [~, brain_outline] = plotBrainGrid([], h2);
                hold on;
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
                [projData.projection_volume(theseTargets)]*100, cmapData(round(round(transp, 2).*100), :), 'filled', ...
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
    %CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_str_to_gpe'], OptionZ)

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

        %load([JF_Scripts_CortexLabPath, 'queryAllenAtlas/Str_snr_medial.mat'])

        for iView = 1 %:3
            h2 = subplot(1, 1, iView);
            if iInjection == size(injectionAreas, 2)
                [~, brain_outline] = plotBrainGrid([], h2);
                hold on;
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
                [projData.projection_volume(theseTargets)]*100, cmapData(round(round(transp, 2).*100), :), 'filled', ...
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
    %CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/vid_str_to_snr'], OptionZ)

    z_rel = max_voxel_z - ml_max / 2;
    max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);

    figure();
    subplot(144)
    colormap(cmapData)
    h = colorbar;
    h.Title.String = ['Normalized', newline, 'projection', newline, 'intensity'];
    makepretty;
    set(gcf, 'Color', 'white')
else

    %% plot data slice by slice
    theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};
    figure(1);
    clf;
    figHandle = figure('Name', regionType);
    figHandle.Color = 'w';
    [tv, av, st, bregma] = ya_loadAllenAtlas(atlasBrainRegLocation); % QQ bregma wrong in brainreg atlas


    structure_alpha = 0.2;


    cl_myPaths;


    regions = {'CP', 'GPe', 'SNr'}; % 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
    warning off;
    regions_id = [672, 1022, 381];

    nChunks = 20;
    nBins = 30;
    clearvars chunks_region

    region_ap_boundaries = [152, 307; ...
        219, 287; ...
        311, 363];
    for iRegion = 1:size(regions, 2)
        chunks_region(iRegion, :) = region_ap_boundaries(iRegion, 1):(region_ap_boundaries(iRegion, 2) - region_ap_boundaries(iRegion, 1)) / nChunks:region_ap_boundaries(iRegion, 2);
        chunks_region_ap(iRegion, :) = (bregma(1) - chunks_region(iRegion, :)) * 25;
    end


    for iRegion = 1:size(regions, 2)

        curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));

        figure(1);

        projection_views = repmat([1, 2], nChunks, 1); %[1,2; 1,3; 2,3];%ML/AP, ML/DV, AP/DV


        % initialize variables
        boundary_projection = cell(3, 1);
        projection_view_bins = cell(3, 1);
        projection_view_lims = nan(3, 2, 2);

        for iChunk = 1:nChunks
            clearvars regionLocation
            % get structure boundaries and plot outline
            region_area = permute(av(round(chunks_region(iRegion, iChunk)):1:round(chunks_region(iRegion, iChunk+1)), ...
                1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
            % AP, DV, ML -> ML, AP, DV

            [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
                = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
            thisChunk_AP = regionLocation(projection_views(1, 1), :);
            thisChunk_DV = regionLocation(projection_views(1, 2), :);

            subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
            boundary_projection{iChunk} = boundary(thisChunk_AP', ...
                thisChunk_DV', 0);

            plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk})./10, ...
                regionLocation(projection_views(1, 2), boundary_projection{iChunk})./10, ...
                'Color', theseColors{iRegion});
            axis equal
            axis square
            axis image
            %makepretty;

            projection_view_lims(iChunk, 1, :) = xlim .* 10;
            projection_view_lims(iChunk, 2, :) = ylim .* 10;
            projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
                (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / nBins: ...
                projection_view_lims(iChunk, 1, 2), ...
                projection_view_lims(iChunk, 2, 1): ...
                (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / nBins: ...
                projection_view_lims(iChunk, 2, 2)};

        end


        prettify_plot('XLimits', 'col');


        curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));

        %figure(figHandle);

        projection_views = repmat([1, 3], nChunks, 1); %[1,2; 1,3; 2,3];%ML/AP, ML/DV, AP/DV


        % initialize variables
        boundary_projection = cell(3, 1);
        projection_view_bins = cell(3, 1);
        projection_view_lims = nan(3, 2, 2);

        for iChunk = 1:nChunks
            clearvars regionLocation
            % get structure boundaries and plot outline
            region_area = permute(av(round(chunks_region(iRegion, iChunk)):1:round(chunks_region(iRegion, iChunk+1)), ...
                1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
            % AP, DV, ML -> ML, AP, DV

            [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
                = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
            regionLocation(2, :) = regionLocation(2, :) + round(chunks_region(iRegion, iChunk)) - 1;
            thisChunk_AP = regionLocation(projection_views(1, 1), :);
            thisChunk_DV = regionLocation(projection_views(1, 2), :);

            subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
            boundary_projection{iChunk} = boundary(thisChunk_AP', ...
                thisChunk_DV', 0);

            plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk})./10, ...
                regionLocation(projection_views(1, 2), boundary_projection{iChunk})./10, ...
                'Color', theseColors{iRegion});
            axis equal
            axis square
            axis image

            projection_view_lims(iChunk, 1, :) = xlim .* 10;
            projection_view_lims(iChunk, 2, :) = ylim .* 10;
            projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
                (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / 15: ...
                projection_view_lims(iChunk, 1, 2), ...
                projection_view_lims(iChunk, 2, 1): ...
                (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / 15: ...
                projection_view_lims(iChunk, 2, 2)};

        end


        % prettify_plot;
        if iRegion == 1
            projData = visProjectionData;
        elseif iRegion == 2
            projData = cp_gpe_ProjectionData;
        elseif iRegion == 3
            projData = cp_snr_ProjectionData;
        end

        % centers: Nx3 matrix where each row defines the [x, y, z] center of a sphere
        % diameters: N-element vector defining the diameter of each sphere
        % imageSize: [width, height, depth] of the resulting 3D image

        % Initialize the 3D image

        theseLocations = [projData.max_voxel_z; projData.max_voxel_y; projData.max_voxel_x] ./ 25; % 2.5 (allen to brain reg) * 10 (slices)

        theseLocationsBregmaAbs = theseLocations';
        bregma_ml_point = bregma(3);

        theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1) > bregma_ml_point), 1) = ...
            bregma_ml_point - ((theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1) > bregma_ml_point), 1)) - bregma_ml_point); % squash right hemisphere on the left


        % Initialize the 3D matrix to store filled bins

        nChunks = length(projection_view_bins); % Assuming each chunk has its X, Y bins defined
        projectionMatrix = {};
        nMatrix = {};
        % voxels are 25um * 10um * 10um, projection volume is in mm^3
        % -> x 2.5*10^6 to convert  radius = projData.projection_volume(iProj) *2.5*10^6; % Influence radius
        % Loop over each chunk
        for iChunk = 1:nChunks
            nBinsX = length(projection_view_bins{iChunk}{1}); % - 1; % Number of bins in X, assuming uniform across chunks
            nBinsY = length(projection_view_bins{iChunk}{2}); % - 1; % Number of bins in Y, assuming uniform across chunks

            projectionMatrix{iChunk} = zeros(nBinsX, nBinsY);
            nMatrix{iRegion}{iChunk} = zeros(nBinsX, nBinsY);
            % Retrieve bin edges for the current chunk
            xEdges = projection_view_bins{iChunk}{1};
            yEdges = projection_view_bins{iChunk}{2};


            % Calculate the size of bins based on the edges
            binSizeX = mean(diff(xEdges));
            binSizeY = mean(diff(yEdges));

            % Given volume in cubic millimeters
            V_um3 = projData.projection_volume ./ projData.volume .* 2.5 .* 1000; % Example volume; replace 100 with your actual volume

            % Calculate the radius in micrometers
            r_um = ((3 * V_um3) / (4 * pi)).^(1 / 3);


            % Iterate through each projection
            % Find unique strings and their indices
            [uniqueStrings, ~, idx] = unique(projData.injection_area);
            % Count the occurrences of each unique string
            occurrences = histcounts(idx, 1:length(uniqueStrings)+1);

            for iProj = 1:size(theseLocationsBregmaAbs, 1)
                thisProj_injectionArea = projData.injection_area{iProj};
                thisProj_injectionAreaN = occurrences(strcmp(thisProj_injectionArea, uniqueStrings));
                % Convert the radius from micrometers to the number of bins it spans
                radiusMicrometers = r_um(iProj); % Adjusted radius in micrometers
                nBinsX_proj = round(radiusMicrometers/binSizeX); % Convert radius to bins in X
                nBinsY_proj = round(radiusMicrometers/binSizeY); % Convert radius to bins in Y

                % Use histcounts2 to find the bin index for the current projection
                [N, ~, ~, binX, binY] = histcounts2(theseLocationsBregmaAbs(iProj, 1), theseLocationsBregmaAbs(iProj, 2), xEdges, yEdges);

                if ~isempty(binX(binX ~= 0)) && ~isempty(binY(binY ~= 0))
                    % Ensure the indices and radius do not exceed the matrix bounds
                    for i = 1:length(binX)
                        if binX(i) > 0 && binY(i) > 0
                            % Calculate the range, ensuring it stays within matrix bounds
                            xRange = max(1, binX(i)-nBinsX_proj):min(nBinsX, binX(i)+nBinsX_proj);
                            yRange = max(1, binY(i)-nBinsY_proj):min(nBinsY, binY(i)+nBinsY_proj);

                            % Increment the intensity within the radius of influence
                            projectionMatrix{iChunk}(xRange, yRange) = ...
                                projectionMatrix{iChunk}(xRange, yRange) + (projData.projection_intensity(iProj) ./ thisProj_injectionAreaN);
                            nMatrix{iRegion}{iChunk}(xRange, yRange) = ...
                                nMatrix{iRegion}{iChunk}(xRange, yRange) + 1;
                        end
                    end
                end
            end
            %projectionMatrix{iChunk} = projectionMatrix{iChunk}./nMatrix{iRegion}{iChunk};
        end

        % for iChunk=1:nChunks
        % subplot(1,nChunks,iChunk)
        % imagesc(projectionMatrix{iChunk}(:,:))
        % end

        % Dimensions for the 3D matrix
        % Initialize parameters
        nChunks = length(projection_view_bins); % Assuming each chunk has its X, Y bins defined
        nBinsX = length(projection_view_bins{1}{1}); % Number of bins in X, assuming uniform across chunks
        nBinsY = length(projection_view_bins{1}{2}); % Number of bins in Y, assuming uniform across chunks

        % normalize each slice
        for iChunk = 1:nChunks
            % Extract the current chunk
            currentChunk = projectionMatrix{iChunk}(:, :);

            % Find the minimum and maximum values in the chunk
            minValue = min(currentChunk(:));
            maxValue = max(currentChunk(:));

            % Avoid division by zero if the chunk is uniform
            if maxValue ~= minValue
                % Normalize the chunk
                normalizedChunk = (currentChunk - minValue) / (maxValue - minValue);
            else
                % If the chunk is uniform, consider it as all zeros or handle as needed
                normalizedChunk = zeros(size(currentChunk));
            end

            % Place the normalized chunk back into the matrix
            projectionMatrix{iChunk}(:, :) = normalizedChunk;
        end

        % % Initialize the 3D matrix to store filled bins
        % projectionMatrix_s = zeros(nBinsX, nBinsY, nChunks);
        %
        % % Iterate through each chunk to fill the matrix based on projection data
        % for iChunk = 1:nChunks
        %     % Retrieve the bin edges for the current chunk
        %     x_binedgesize=mean(diff((projection_view_bins{iChunk}{1})));
        %     xEdges = min(projection_view_bins{iChunk}{1})-(50*x_binedgesize):x_binedgesize:max(projection_view_bins{iChunk}{1})+(50*x_binedgesize);
        %      y_binedgesize=mean(diff((projection_view_bins{iChunk}{2})));
        %     yEdges = min(projection_view_bins{iChunk}{2})-(50*y_binedgesize):y_binedgesize:max(projection_view_bins{iChunk}{2})+(50*y_binedgesize);
        %
        %     % Use histcounts2 to categorize each projection into bins based on its X and Y coordinates
        %     [N, ~, ~, binX, binY] = histcounts2(theseLocationsBregmaAbs(:, 1), theseLocationsBregmaAbs(:, 2), xEdges, yEdges);
        %
        %     % Fill the corresponding bins in the projectionMatrix with 1s
        %     % N contains the counts of projections per bin, binX and binY are the bin assignments for each projection
        %     for i = 1:length(binX)
        %         if binX(i) > 0 && binY(i) > 0
        %             projectionMatrix_s(binX(i), binY(i), iChunk) = 1;
        %         end
        %     end
        % % end
        %
        %         figure();
        %         for iChunk = 1:nChunksfigure(figHandle);
        %             subplot(1, nChunks, iChunk)
        %             imagesc(projectionMatrix_s(:, :, iChunk))
        %         end


        % smooth
        dataSmoothed = cell(nChunks, 1); % Initialize cell array for smoothed data

        % Define Gaussian kernel for smoothing
        gaussianSize = 8; % Size of the Gaussian filter
        gaussianSigma = 4; % Standard deviation of the Gaussian filter
        G = fspecial('gaussian', gaussianSize, gaussianSigma);

        % Apply smoothing to each slice
        for i = 1:nChunks
            currentSlice = projectionMatrix{i};
            if ~isempty(currentSlice) % Ensure the slice is not empty
                % Apply Gaussian smoothing
                smoothedSlice = imfilter(currentSlice, G, 'same', 'replicate');
                dataSmoothed{i} = smoothedSlice;
            else
                dataSmoothed{i} = currentSlice; % In case of an empty slice, just carry it over
            end
        end
        % smooth across slices
        nChunks = numel(dataSmoothed);


        % Here, we focus on smoothing across the Z dimension (between slices)

        % Parameters for 1D Gaussian smoothing across slices
        zGaussianSize = 2; % Smaller size as we're only smoothing between a few slices
        zGaussianSigma = 1; % Adjust based on the desired level of cross-slice smoothing
        zG = fspecial('gaussian', [zGaussianSize, 1], zGaussianSigma); % 1D Gaussian

        % Initialize an array to store the smoothed data (using the largest slice as a reference)
        % For simplicity, let's assume we've interpolated all slices to a common size (this step is not shown here)
        % maxRows and maxCols represent these common dimensions
        [maxRows, maxCols] = size(dataSmoothed{1}); % Assuming dataSmoothed{1} is your reference size
        tempVolume = zeros(maxRows, maxCols, nChunks);

        % Populate the temporary 3D volume with data from smoothed sprojection_view_binslices
        for i = 1:nChunks
            tempSlice = dataSmoothed{i};
            % You might need to interpolate tempSlice to [maxRows, maxCols] here if sizes vary
            tempVolume(:, :, i) = tempSlice;
        end

        % Apply 1D Gaussian smoothing across the Z dimension (between
        % slices) -> disbaled for now
        % for x = 1:maxRows
        %     for y = 1:maxCols
        %         pixelSeries = squeeze(tempVolume(x, y, :)); % Extract the pixel series across slices
        %         smoothedPixelSeries = imfilter(pixelSeries, zG, 'same', 'replicate'); % Apply 1D smoothing
        %         tempVolume(x, y, :) = smoothedPixelSeries; % Update the volume with smoothed values
        %     end
        % end

        % Extract the smoothed slices back into the cell array format
        dataCrossSmoothed = cell(nChunks, 1);
        for i = 1:nChunks
            dataCrossSmoothed{i} = tempVolume(:, :, i);
        end


        % figure();
        % for iChunk = 1:nChunks
        %     subplot(1, nChunks, iChunk)
        %     imagesc(dataCrossSmoothed{iChunk}(:, :))
        % end
for iChunk=1:nChunks 
    maxValue_chunks(iChunk) = nanmax(nanmax(dataCrossSmoothed{iChunk}));
end
maxValue = max(maxValue_chunks);
        %% plot average increase for each bin

        for iChunk = 1:nChunks
            %maxValue = nanmax(nanmax(dataCrossSmoothed{iChunk}));


            thisCmap_limits = [-maxValue, maxValue];

            binnedArrayPixelSmooth = dataCrossSmoothed{iChunk};
            % remove any data points outside of the ROI
            clearvars regionLocation
            % get structure boundaries and plot outline
            region_area = permute(av(round(chunks_region(iRegion, iChunk)):1:round(chunks_region(iRegion, iChunk+1)), ...
                1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
            % AP, DV, ML -> ML, AP, DV

            [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
                = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
            isIN = nan(size(binnedArrayPixelSmooth, 1), size(binnedArrayPixelSmooth, 2));
            for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
                for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
                    isIN(iPixelX, iPixelY) = inpolygon(projection_view_bins{iChunk}{1}(iPixelX), ...
                        projection_view_bins{iChunk}{2}(iPixelY), ...
                        regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
                        regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}));
                end
            end

            figure(figHandle);
            subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
            %if pcells
            binnedArrayPixelSmooth(isIN == 0) = NaN;
            %binnedArrayPixelSmooth(binnedArrayPixelSmooth == 0) = 0.0001; %convert 0s to very small values - just for plotting purposes
            %else
            %    binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
            %end
            ax = gca;
            ax.YColor = 'w'; % Red
            ax.XColor = 'w'; % Red
            im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
                binnedArrayPixelSmooth');
            set(im, 'AlphaData', ~isnan(get(im, 'CData')));

            set(gca, 'color', [0.5, 0.5, 0.5]);

            originalColormap = brewermap([], '*RdBu');

            colormap(originalColormap);


            try
                caxis(thisCmap_limits)
            catch
                caxis([-1, 1])
            end
            hold on;

            %colorbar
            clearvars binnedArrayPixel
            hold on;
            boundary_projection{iChunk} = boundary(regionLocation(projection_views(iChunk, 1), :)', ...
                regionLocation(projection_views(iChunk, 2), :)', 0);
            plot(regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
                regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}), ...
                'Color', theseColors{iRegion}, 'LineWidth', 2);

            axis equal
            axis square
            axis image
            ax.XLabel.Color = [0, 0, 0];
            ax.YLabel.Color = [0, 0, 0];
            nColors = numel(ax.YTickLabel);
            cm = [0, 0, 0];
            for i = 1:nColors
                ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
            end

            nColors = numel(ax.XTickLabel);
            cm = [0, 0, 0];
            for i = 1:nColors
                ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
            end

            %prettify_plot;
            clearvars isIN
            try
                caxis(thisCmap_limits)
            catch
                caxis([-1, 1])
            end
            set(gca, 'color', [0.5, 0.5, 0.5]);
            xlim([projection_view_bins{iChunk}{1}(1), projection_view_bins{iChunk}{1}(end)])
            ylim([projection_view_bins{iChunk}{2}(1), projection_view_bins{iChunk}{2}(end)])

            axis off;

            % Alternatively, to remove only specific components:
            set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
            xlabel(''); % Remove x-axis label
            ylabel(''); % Remove y-axis label
            thisAp_slice = round(nanmean(chunks_region_ap(iRegion, iChunk:iChunk+1))./1000, 2);
            coordinates = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
            if iChunk == 1
                
                title(['Bregma: ', num2str(thisAp_slice)]); % Remove title
            elseif iChunk == nChunks
                title([num2str(thisAp_slice), 'mm']); % Remove title
            else
                title([num2str(thisAp_slice)]); % Remove title
            end
            % title(''); % Remove title


        end
        % set x and ylims
        xlims_region = nan(nChunks, 2);
        ylims_region = nan(nChunks, 2);
        for iChunk = 1:nChunks
            subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
            xlims_region(iChunk, :) = xlim;
            ylims_region(iChunk, :) = ylim;
        end

        diff_xlims_region = diff(xlims_region');
        diff_ylims_region = diff(ylims_region');
        for iChunk = 1:nChunks
            subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
            xlims_here = (max(diff_xlims_region) - diff_xlims_region(iChunk)) ./ 2;
            xlim([xlims_region(iChunk, 1) - xlims_here, xlims_region(iChunk, 2) + xlims_here])

            ylims_here = (max(diff_ylims_region) - diff_ylims_region(iChunk)) ./ 2;
            ylim([ylims_region(iChunk, 1) - ylims_here, ylims_region(iChunk, 2) + ylims_here])
            if iChunk == 1
                axis_length_mm = 1;
                one_pixel_x = (diff(projection_view_bins{iChunk}{1}(2:3)));
                one_pixel_x_um = one_pixel_x * 25;
                one_pixel_y = (diff(projection_view_bins{iChunk}{2}(2:3)));
                one_pixel_y_um = one_pixel_y * 25;
                axis_length_atlas_units_x = (axis_length_mm * 1000) / (one_pixel_x_um);
                axis_length_atlas_units_y = (axis_length_mm * 1000) / (one_pixel_y_um);
                prettify_addScaleBars(axis_length_atlas_units_x, axis_length_atlas_units_y, ...
                    [num2str(axis_length_mm), 'mm'], [num2str(axis_length_mm), 'mm'], 'topLeft', '', '')
            end

        end


        keep injectionAreas plotInjections regionType tv figHandle regions thisCmap_limits st av structure_alpha chunks_region_ap theseColors iRegion bregma nChunks chunks_region nBins cp_gpe_ProjectionData cp_snr_ProjectionData visProjectionData
    end
    %prettify_plot;
end
if plotInjections
%% plot injection sites in 3D
cl_myPaths;
cl_plottingSettings;
tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean


injectionTable = readtable([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/DATA/allenAtlasProjection_info.csv']);

% injection coordinates
combinedStr = strcat(injectionTable.injection_coordinates{:});
combinedStr = strrep(combinedStr, ']', ', ');
cleanStr = erase(combinedStr, {'[', ']'});
cleanStr(end-1:end) = '';
numStrs = strsplit(cleanStr, ',');
nums = str2double(numStrs);
injection_coordinatesMatrix = reshape(nums, 3, []).'; % Adjust dimensions as necessary

% Convert cell array to char array for vectorized operations
hexMatrix = char(injectionTable .structure_color);
rHex = hexMatrix(:, 1:2);
gHex = hexMatrix(:, 3:4);
bHex = hexMatrix(:, 5:6);
r = hex2dec(rHex);
g = hex2dec(gHex);
b = hex2dec(bHex);
rgbMatrix = [r, g, b];


theseInjections = ismember(injectionTable.structure_abbrev, injectionAreas(1:end-1));

curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = regionColors{1};

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;


[h1, patchHandle] = ya_plotBrainSurface(allenAtlasPath, 'w');


hold on;
structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);

iView = 1;
axis equal
[ap_max, dv_max, ml_max] = size(tv);

hold on;


max_voxel_x = [injection_coordinatesMatrix(theseInjections, 1)] / 10;
max_voxel_y = [injection_coordinatesMatrix(theseInjections, 2)] / 10;
max_voxel_z = [injection_coordinatesMatrix(theseInjections, 3)] / 10;

structure_color = rgbMatrix(theseInjections, :) ./ 255;

scatter3(max_voxel_x, max_voxel_z, max_voxel_y, ...
    [sqrt(injectionTable.injection_volume(theseInjections))]*50, structure_color, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', ...
    1)
if iView ~= 1
    set(h1, 'Ydir', 'reverse')
end

set(gcf, 'Color', 'w')

%% gpe
expIDs_medialStriatum = [ ...
    124059700, ...
    159223001, ... %anterior
    146553266, ... %keep? anterior
    112458831, ... %keep?
    100142580, ... %keep??
    301620241, ... %keep?
    278434443, ...
    127140981, ... %keep??
    308395312, ... %keep??
    159941339, ...
    307910595, ... %keep??
    112307754, ... %keep?
    293366741, ...
    293366035, ...
    160540013, ...
    301180385, ...
    112458831, ...
    127762867, ... %posterior
    117317884, ... %posterior
    120762196]; %posterior

expIDs_centralStriatum = [307910595, ...
    272697238, ...
    264095536, ...
    112307754, ...
    120570964, ...
    160537796, ...
    127140981, ...
    157911832];


expIDs_medialPFC = [183618139, ...
    112952510, ...
    177459319, ...
    184169615, ...
    166460484, ...
    166055636, ...
    182814071, ...
    292172846, ...
    298000163];

if strcmp(regionType, 'VIS')
    theseInjections = ismember(injectionTable.id, expIDs_medialStriatum);
elseif strcmp(regionType, 'mPFC')
    theseInjections = ismember(injectionTable.id, expIDs_centralStriatum);
end


curr_plot_structure_idx = find(contains(st.name, 'pallidus ex'));
slice_spacing = 10;
plot_structure_color = regionColors{2};

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.5;


[h1, patchHandle] = ya_plotBrainSurface(allenAtlasPath, 'w');


hold on;
structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);

iView = 1;
axis equal
[ap_max, dv_max, ml_max] = size(tv);

hold on;


max_voxel_x = [injection_coordinatesMatrix(theseInjections, 1)] / 10;
max_voxel_y = [injection_coordinatesMatrix(theseInjections, 2)] / 10;
max_voxel_z = [injection_coordinatesMatrix(theseInjections, 3)] / 10;

structure_color = rgbMatrix(theseInjections, :) ./ 255;

scatter3(max_voxel_x, max_voxel_z, max_voxel_y, ...
    [sqrt(injectionTable.injection_volume(theseInjections))]*50, regionColors{1}, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
    1)
if iView ~= 1
    set(h1, 'Ydir', 'reverse')
end

set(gcf, 'Color', 'w')
%% snr


curr_plot_structure_idx = find(contains(st.name, 'ubstantia nigra reti'));
slice_spacing = 10;
plot_structure_color = regionColors{3};

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.5;


[h1, patchHandle] = ya_plotBrainSurface(allenAtlasPath, 'w');


hold on;
structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);

iView = 1;
axis equal
[ap_max, dv_max, ml_max] = size(tv);

hold on;


max_voxel_x = [injection_coordinatesMatrix(theseInjections, 1)] / 10;
max_voxel_y = [injection_coordinatesMatrix(theseInjections, 2)] / 10;
max_voxel_z = [injection_coordinatesMatrix(theseInjections, 3)] / 10;

structure_color = rgbMatrix(theseInjections, :) ./ 255;

scatter3(max_voxel_x, max_voxel_z, max_voxel_y, ...
    [sqrt(injectionTable.injection_volume(theseInjections))]*50, regionColors{1}, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
    1)
if iView ~= 1
    set(h1, 'Ydir', 'reverse')
end

set(gcf, 'Color', 'w')
end
end