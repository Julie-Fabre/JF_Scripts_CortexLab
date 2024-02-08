
%% cl_locationPSTH
passive_data = load('/home/julie/Dropbox/MATLAB/naive_data1.mat');
regions = {'CP', 'GPe', 'SNr'};

keep passive_data regions
figure(1);clf;
figure(2);clf;

cl_myPaths;
regionResolution = [1, 1, 1, 1, 1, 1, 1];
% ARA_levels = [47,66; ...% CP
%     58, 68;...%GPe
%     81, 90;...%SNr
%     NaN, NaN;...%GPi
%     NaN,NaN;...%STN
%     81, 90;...%SNr
%     NaN, NaN;...%VTA
%     NaN,NaN];...%SNc;
passive_data.psth_average = squeeze(nanmean(passive_data.psth{1}(:, 3, :, :), 3));
zscore_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
    nanstd(passive_data.psth_average(:, 1:200), [], 2);
dFR_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
    nanmean(passive_data.psth_average(:, 1:200), 2);

thisCmap_limits = [-300, 300];
theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

if ~exist('st', 'var')
    [tv, av, st, bregma] = ya_loadAllenAtlas(atlasBrainRegLocation); % QQ bregma wrong in brainreg atlas
end

structure_alpha = 0.2;

%% get each region's limits and define chunks
nChunks = 4;
clearvars chunks_region
% for iRegion = 1:size(regions, 2)
%     clearvars regionLocation
%     curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));
% 
%     region_area = permute(av(1:1:end, ...
%         1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
%     % AP, DV, ML -> ML, AP, DV
% 
%     [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
%         = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
% 
%     region_ap_boundaries(iRegion, :) = [min(regionLocation(:, 2)), max(regionLocation(:, 2))];
% 
%     chunks_region(iRegion,:) = region_ap_boundaries(iRegion, 1):(region_ap_boundaries(iRegion, 2)-region_ap_boundaries(iRegion, 1))/nChunks:region_ap_boundaries(iRegion, 2);
% end
region_ap_boundaries = [152, 307;...
    219,287;...
    311,363];

for iRegion = 1:size(regions, 2)
    chunks_region(iRegion,:) =  region_ap_boundaries(iRegion, 1):(region_ap_boundaries(iRegion, 2)-region_ap_boundaries(iRegion, 1))/nChunks:region_ap_boundaries(iRegion, 2);
end

%% Visualize AP chunks in ML x AP projection
for iRegion = 1:size(regions, 2)

    curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));

    figure(1);

    projection_views = repmat([1,2],nChunks, 1);%[1,2; 1,3; 2,3];%ML/AP, ML/DV, AP/DV


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

        subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (4 - iChunk)))
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
            (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / 15: ...
            projection_view_lims(iChunk, 1, 2), ...
            projection_view_lims(iChunk, 2, 1): ...
            (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / 15: ...
            projection_view_lims(iChunk, 2, 2)};

    end

    theseLocations = passive_data.unit_coords;
    theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
        theseLocations(:, 1), ...
        theseLocations(:, 2)]; %AP, DV, ML -> ML, AP, DV
    
    prettify_plot('XLimits', 'col');
end


%% AP chunks in ML x DV projection 
for iRegion = 1:size(regions, 2)

    curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));

    figure(2);

    projection_views = repmat([1,3],nChunks, 1);%[1,2; 1,3; 2,3];%ML/AP, ML/DV, AP/DV


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

        subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (4 - iChunk)))
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

     theseLocations = passive_data.unit_coords;
    theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
        theseLocations(:, 1), ...
        theseLocations(:, 2)];% go from AP, DV, ML to ML, AP, DV (like loaded Atlas) 

    bregma_ml_point = bregma(1) / 2.5; %2.5 is difference in scaling between 
    % brainreg (25 um resolution) and allen (10um resolution, where this bregma value comes from)

    theseLocationsBregmaAbs(:,1) = bregma_ml_point - abs(theseLocationsBregmaAbs(:,1) - bregma_ml_point); % squash right hemisphere on the left

    %% plot average increase for each bin

    for iChunk = 1:nChunks
        [N, Xedges, Yedges, binX, binY] = histcounts2(theseLocationsBregmaAbs(:, projection_views(iChunk, 1)), ...
            theseLocationsBregmaAbs(:, projection_views(iChunk, 2)), projection_view_bins{iChunk}{1}, ...
            projection_view_bins{iChunk}{2}); %par rapport a bregma!

        binnedArrayPixel = nan(size(Xedges, 2), size(Yedges, 2)); % initialize

        for iBinX = 1:size(Xedges, 2)
            for iBinY = 1:size(Yedges, 2)
                theseNeurons = binX == iBinX & binY == iBinY & passive_data.unit_area == iRegion &...
                (passive_data.unitType' ==1 | passive_data.unitType' ==2);
                binnedArrayTot = [];
                if sum(theseNeurons) > 0
                    %remove any infs QQ i need to deal with this!
                    mean_2d = nanmean(abs(dFR_psth(theseNeurons, :)), 2);
                    binnedArrayPixel(iBinX, iBinY) = nanmean(mean_2d(~isinf(mean_2d)));

                end

            end
        end

        %NaN-out 0 values (= no data) and Inf values if present (QQ need to check why Inf values sometimes crop up)
        binnedArrayPixel(binnedArrayPixel == Inf) = NaN;
        binnedArrayPixel(binnedArrayPixel == 0) = NaN;

        % smooth data
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 4, 4);
        %binnedArrayPixelSmooth = binnedArrayPixel;

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

        figure(2);
        subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (iChunk - 1)))

        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
        ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            binnedArrayPixelSmooth'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));

        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        caxis(thisCmap_limits)
        hold on;

        %colorbar
        clearvars binnedArrayPixel
        hold on;
        boundary_projection{iChunk} = boundary(regionLocation(projection_views(iChunk, 1), :)', ...
            regionLocation(projection_views(iChunk, 2), :)', 0);
        plot(regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
            regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}), ...
            'Color', theseColors{iRegion});

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
        caxis(thisCmap_limits)
        set(gca, 'color', [0.5, 0.5, 0.5]);
        xlim([projection_view_bins{iChunk}{1}(1), projection_view_bins{iChunk}{1}(end)])
        ylim([projection_view_bins{iChunk}{2}(1), projection_view_bins{iChunk}{2}(end)])
    end
    keep passive_data regions thisCmap_limits st av regionResolution structure_alpha theseColors dFR_psth iRegion bregma nChunks chunks_region
end

