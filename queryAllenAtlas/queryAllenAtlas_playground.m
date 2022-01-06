%docs: https://allensdk.readthedocs.io/en/latest/unionizes.html
%% info
myPaths;
allen_atlas_path = [allenAtlasPath 'allenCCF'];
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;

curr_plot_structure = st.id(find(contains(st.name, 'audoputamen')));

%% get all VIS Ctx experiments
injectionAreas = {'VISp', 'VISl', 'VISpl', 'VISpm', 'VISal', 'VISam'};
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

for iInjection = 1:size(injectionAreas, 2)

    if ~exist([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'], 'file')
        expIDs = findAllenExperiments('injection', injectionAreas{iInjection}, 'primary', true);
        proj = [];
        for iExpID = 1:size(expIDs, 2)
            proj_temp = getProjectionDataFromExperiment(expIDs(iExpID));
            proj = [proj, proj_temp{1, 1}];
        end
        save([dropboxPath, 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'], 'proj')
    else
        load([dropboxPath 'MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/', injectionAreas{iInjection}, '.mat'])
    end

    theseTargets = ismember([proj.structure_id], double(curr_plot_structure)) & [proj.max_voxel_density] > 0;

    visProjectionData.max_voxel_x = [visProjectionData.max_voxel_x, [proj(theseTargets).max_voxel_x]];
    visProjectionData.max_voxel_y = [visProjectionData.max_voxel_y, [proj(theseTargets).max_voxel_y]];
    visProjectionData.max_voxel_z = [visProjectionData.max_voxel_z, [proj(theseTargets).max_voxel_z]];
    visProjectionData.normalized_projection_volume = [visProjectionData.normalized_projection_volume, [proj(theseTargets).normalized_projection_volume]];
    visProjectionData.projection_intensity = [visProjectionData.projection_intensity, [proj(theseTargets).projection_intensity]];
    areaName = cell(1, size(find(theseTargets), 2));
    areaName(:) = {injectionAreas{iInjection}};
    visProjectionData.injection_area = [visProjectionData.injection_area, areaName];

    clearvars areaName proj


end
save([dropboxPath, '/MATLAB/onPaths/JF_Scripts_CortexLab/queryAllenAtlas/VIS_Str.mat'], 'visProjectionData')

%% plot figure 3D
figure();
for iInjection = 1:size(injectionAreas, 2)

    load([JF_Scripts_CortexLabPath 'queryAllenAtlas/VIS_Str.mat'])

    for iView = 1:3
        h1 = subplot(1, 4, iView);
        %[~, brain_outline] = plotBrainGrid([], ss);
        hold on;
        axis equal
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        if collapseHemispheres
            ylim([-10, ml_max / 2 + 10])
        else
            ylim([-10, ml_max + 10])
        end

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
%         structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
%             'Faces', structure_3d.faces, ...
%             'FaceColor', 'w', 'EdgeColor', 'None', 'FaceAlpha', 1);
        %         minV = min(min(structure_patch.Vertices(:,:)));
        %         maxV = max(max(structure_patch.Vertices(:,:)));
        %         xlim([])


        theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        transp = [projData.projection_intensity(theseTargets)] ./ (max([projData.projection_intensity(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt==0)=1;
        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [projData.normalized_projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
            1)
        if iView ~= 1
            set(h1, 'Ydir', 'reverse')
        end
set(h1, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
xlim([300 800])
set(gca,'Visible','off')
makepretty;
title(injectionAreas{iInjection})
    end
end
subplot(144)
colormap(cmapData)
h=colorbar;
h.Title.String = ['Normalized' newline 'projection' newline 'intensity'];
makepretty;
set(gcf, 'Color', 'white')

%% plot figure 2D
iRegion=1;
figure();
for iInjection = 1:size(injectionAreas, 2)

    load([JF_Scripts_CortexLabPath 'queryAllenAtlas/VIS_Str.mat'])

       curr_plot_structure = find(strcmp(st.acronym, regions{iRegion}));
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);

        ii = permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
        ffind = find(ii);
        [r, c, v] = ind2sub(size(ii), [ffind]);

        figure(1);

        XYZ = {r, c; r, v; c, v};
        XYZnormTerm = [-(bregma(3) / 100),  - (bregma(1) / 100); -(bregma(3) / 100), 0 ; - (bregma(1) / 100), 0]; 
        XYZnormFactor = [1,-1; 1, -1; -1, -1]; 
        xyzInd = [3,1;2,1;2,3];
    for iView = 1:3
        iXYZ=iView;
        subplot(1,4, iView)
            cla;
            bb = boundary(XYZ{iXYZ,1},XYZ{iXYZ,2}, 0);
            xx = XYZnormFactor(iXYZ,1)*(XYZ{iXYZ,1}(bb)./10+XYZnormTerm(iXYZ,1));
            yy = XYZnormFactor(iXYZ,2)*(XYZ{iXYZ,2}(bb)./10+XYZnormTerm(iXYZ,2));
            
            thisCentroid = [mean(xx),mean(yy)];
            
            
            plot(xx, yy, 'Color', theseColors{iRegion});
           
            %hold on;
            axis equal
      

            xyXLim = xlim;
            xyYLim = ylim;
            makepretty;
             yy_new=yy;
            xx_new=xx;
              xyzCountBins(iXYZ,:) = {xyXLim(1) * 10 :(xyXLim(2) - xyXLim(1)) * 10 / 15:xyXLim(2) * 10, ...
                xyYLim(1) * 10 :(xyYLim(2) - xyYLim(1)) * 10 / 15:xyYLim(2) * 10};
            newXXYY{iXYZ,:,:}=[xx_new, yy_new];
            plot(newXXYY{iXYZ,:,:}(:,1), newXXYY{iXYZ,:,:}(:,2), 'Color', theseColors{iRegion},'LineWidth',5);
            xlabel('ML (mm)')
axis equal
        %zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        %         if iInjection == 1
        %             xlabel(viewLegends{iView, 1})
        %             ylabel(viewLegends{iView, 2})
        %         end

        % xticks([0 5 10])
        % xticklabels({'x = 0','x = 5','x = 10'})
        projData = visProjectionData;
        [ap_max, dv_max, ml_max] = size(tv);
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
bregma = [540,0,570];
        theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        transp = [projData.projection_intensity(theseTargets)] ./ (max([projData.projection_intensity(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt==0)=1;
%         if iXYZ == 3 
%              scatter(max_voxel_x(theseTargets)./100, max_voxel_y(theseTargets)./100, ...
%             [projData.normalized_projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
%             'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
%             1)
%         elseif iXYZ == 2
%              scatter( max_voxel_y(theseTargets)./100, max_voxel_z(theseTargets)./100, ...
%             [projData.normalized_projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
%             'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
%             1)
%         elseif iXYZ == 1
%              scatter(max_voxel_x(theseTargets)./100 - (bregma(1)./100) -3, (bregma(3)./100) - max_voxel_z(theseTargets)./100 -3,  ...
%             [projData.normalized_projection_volume(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
%             'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
%             1)
%         end


set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
%xlim([300 800])
set(gca,'Visible','off')
makepretty;
title(injectionAreas{iInjection})
    end
end
subplot(144)
colormap(cmapData)

makepretty;
set(gcf, 'Color', 'white')

%% DMS to GPe

curr_plot_structure = st.id(find(contains(st.name, 'obus pallidus ext')));
result = getProjectionDataFromExperiment(124059700);
projData = result{1,1};
figure();
theseTargets = ismember([projData.structure_id], double(curr_plot_structure)) & [projData.max_voxel_density] > 0;
for iView = 1:3
    h1 = subplot(1, 4, iView);
    view(viewVals(iView, :))
     hold on;
        axis equal
        view(viewVals(iView, :))
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        if collapseHemispheres
            ylim([-10, ml_max / 2 + 10])
        else
            ylim([-10, ml_max + 10])
        end

        zlim([-10, dv_max + 10])
        %in tru coordinates
        hold on;
        
max_voxel_x = [projData.max_voxel_x] / 10;
            max_voxel_y = [projData.max_voxel_y] / 10;
            max_voxel_z = [projData.max_voxel_z] / 10;
 max_voxel_z = [projData.max_voxel_z] / 10;
            z_rel = max_voxel_z - ml_max / 2;
            max_voxel_z(z_rel > 0) = ml_max / 2 - z_rel(z_rel > 0);

            max_voxel_y = [projData.max_voxel_y] / 10;
             % theseTargets = strcmp(projData.injection_area, injectionAreas(iInjection));
        hold on;
        qq=[projData.projection_intensity];
        transp = [qq(theseTargets)] ./ (max([qq(theseTargets)]));
        tt = round(round(transp, 2).*100);
        tt(tt==0)=1;
        ww=[projData.normalized_projection_volume];
        scatter3(max_voxel_x(theseTargets), max_voxel_z(theseTargets), max_voxel_y(theseTargets), ...
            [ww(theseTargets)]*10, cmapData(round(round(transp, 2).*100), :), 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', ...
            1)
               set(h1, 'Ydir', 'reverse')
        
%set(h1, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
xlim([300 800])
%set(gca,'Visible','off')
makepretty;
title(injectionAreas{iInjection})

end
%% PS to GPe 