% to do 
% add probe depths, types, and mutli color if several regions of interest
% plot only probe in the region of interest? 

animalsType = {'MultiImgWorld'};
regionsNames = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
regions = {'DMS', 'STN', 'GPe', 'SNr', 'GPi'};
recordingInfo = readtable('C:\Users\Julie\Dropbox\Analysis\Recordings - Sheet1.csv');

%add probe types, depths
for iType = 1:size(animalsType, 2)
    theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    theseColors = {rgb('DeepSkyBlue');  rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black');rgb('Brown')};

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
    for iRegion = 1:5
        curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
         
        if strcmp(regionsNames{iRegion}, 'STR') == 0
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
        end
        %plot probe tracks
        theseProbes = contains(recordingInfo.Location, regions{iRegion});
        %get animal and probe, load track
        theseAnimals = recordingInfo.Mouse(theseTypes);

        for iAnimal = [1,2,4]%1:size(theseAnimals, 1)
            %iAnimal = iAnimal + 1;
            load(['D:\', theseAnimals{iAnimal}, '\slices\probe_ccf.mat'])
            qqqq=unique(theseAnimals);
            thisAnimal = strcmp(recordingInfo.Mouse, qqqq{iAnimal});
            ttP = (recordingInfo.HistologyProbe( thisAnimal & theseProbes));
            %thisProbe = recordingInfo.HistologyProbe(find(thisAnimal & theseProbes));
            for iProbe = 1:size(probe_ccf,1)
                %iProbe = iProbe +1 
                thisthisProbe = iProbe;
                if (any(probe_ccf(iProbe).trajectory_areas == curr_plot_structure) ) || (iAnimal == 2 && ismember(iProbe,ttP))
                r0 = mean(probe_ccf(thisthisProbe).points, 1);
                xyz = bsxfun(@minus, probe_ccf(thisthisProbe).points, r0);
                [~, ~, V] = svd(xyz, 0);
                histology_probe_direction = V(:, 1);

                probe_eval_points = [-1000, 1000];
                probe_line_endpoints = bsxfun(@plus, bsxfun(@times, probe_eval_points', histology_probe_direction'), r0);

                % Philip's GUI: not saved in native CCF order?
                % plot3(probe_ccf(iProbe).points(:,3),probe_ccf(iProbe).points(:,1),probe_ccf(iProbe).points(:,2),'.b','MarkerSize',20);
                % line(P(:,3),P(:,1),P(:,2),'color','k','linewidth',2)

                % % Mine: saved in native CCF order [AP,DV,ML]
                plot3(probe_ccf(thisthisProbe).points(:, 1), probe_ccf(thisthisProbe).points(:, 3), probe_ccf(thisthisProbe).points(:, 2),...
                    '.', 'MarkerSize', 20,'Color', theseColors{iRegion, :});
                line(probe_line_endpoints(:, 1), probe_line_endpoints(:, 3), probe_line_endpoints(:, 2), 'color', ...
                    theseColors{iRegion, :}, 'linewidth', 2)

                % Place the probe on the histology best-fit axis
                [ap_max, dv_max, ml_max] = size(tv);

                probe_ref_top = probe_line_endpoints(1, [1, 3, 2]);
                probe_ref_bottom = probe_line_endpoints(2, [1, 3, 2]);
                probe_ref_vector = [probe_ref_top; probe_ref_bottom]';
                end
%                 probe_vector = [probe_ref_vector(:, 1), diff(probe_ref_vector, [], 2) ./ ...
%                     norm(diff(probe_ref_vector, [], 2)) * probe_length + probe_ref_vector(:, 1)];
            end
        end

    end
end
