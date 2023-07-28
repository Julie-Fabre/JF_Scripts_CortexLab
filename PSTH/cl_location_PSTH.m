
%% cl_locationPSTH

keep passive_data regions
close all;
figure('Color', 'white')
cl_myPaths;
regionResolution = [1, 1, 1, 1, 1, 1, 1];

zscore_psth = (passive_data.psth - nanmean(passive_data.psth(:, 1:50), 2)) ./ ...
    nanstd(passive_data.psth(:, 1:50), [], 2);
dFR_psth = (passive_data.psth - nanmean(passive_data.psth(:, 1:50), 2)) ./ ...
    nanmean(passive_data.psth(:, 1:50), 2);

thisCmap_limits = [-150, 150];
theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

if ~exist('st', 'var')
    [tv, av, st, bregma] = bd_loadAllenAtlas(atlasBrainRegLocation); % QQ bregma wrong in brainreg atlas
end

structure_alpha = 0.2;

%% X-Y, X-Z, Y-Z plots for each region.

for iRegion = 1:size(regions, 2)

    curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));
%     structure_3d = isosurface(permute(av(1:regionSpacing(iRegion):end, ...
%         1:regionSpacing(iRegion):end, 1:regionSpacing(iRegion):end) == curr_plot_structure, [3, 1, 2]), 0);

    % get structure area 
    region_area = permute(av(1:1:end, ...
        1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
    [regionLocation(1,:), regionLocation(2,:), regionLocation(3,:)]...
        = ind2sub(size(region_area), find(region_area)); %ML, AP, DV

    figure(1);
    
    projection_views = [2,1; 2,3; 1,3];

    % initialize variables 
    boundary_projection = cell(3,1);
    projection_view_bins = cell(3,1);
    projection_view_lims = nan(3,2,2);

    for iProjection = 1:3 
        % get structure boundaries and plot outline
        subplot(3, size(regions, 2), iRegion + (size(regions,2) * (iProjection-1)))
        boundary_projection{iProjection} = boundary(regionLocation(projection_views(iProjection,1),:)',...
            regionLocation(projection_views(iProjection,2),:)', 0);
        plot(regionLocation(projection_views(iProjection,1),boundary_projection{iProjection})./10,...
            regionLocation(projection_views(iProjection,2),boundary_projection{iProjection})./10,...
            'Color', theseColors{iRegion});
        axis equal
        axis square
        axis image
        makepretty;

        projection_view_lims(iProjection,1,:) = xlim .* 10;
        projection_view_lims(iProjection,2,:) = ylim .* 10;
        projection_view_bins{iProjection} = {projection_view_lims(iProjection,1,1):...
            (projection_view_lims(iProjection,1,2) - projection_view_lims(iProjection,1,1)) / 15:...
            projection_view_lims(iProjection,1,2),...
            projection_view_lims(iProjection,2,1):...
            (projection_view_lims(iProjection,2,2) - projection_view_lims(iProjection,2,1)) / 15:...
            projection_view_lims(iProjection,2,2)};
        
    end

    theseLocations = passive_data.unit_coords;
    theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
        theseLocations(:, 1), ...
        theseLocations(:, 2)];%AP, DV, ML -> ML, AP, DV

    %% plot average increase for each bin 

    for iProjection = 1:3
        [N, Xedges, Yedges, binX, binY] = histcounts2(theseLocationsBregmaAbs(:, projection_views(iProjection,1)), ...
            theseLocationsBregmaAbs(:, projection_views(iProjection,2)), projection_view_bins{iProjection}{1},...
            projection_view_bins{iProjection}{2}); %par rapport a bregma!

        binnedArrayPixel = nan(size(Xedges, 2), size(Yedges, 2)); % initialize 

        for iBinX = 1:size(Xedges, 2)
            for iBinY = 1:size(Yedges, 2)
                theseNeurons = binX == iBinX & binY == iBinY & passive_data.unit_area == iRegion;% &...
                    %(passive_data.unitType' ==1 | passive_data.unitType' ==2);
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
        %if iRegion==1
            binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 1, 1);
        %else
       %     binnedArrayPixelSmooth = binnedArrayPixel;
        %end

       % remove any data points outside of the ROI
       isIN = nan(size(binnedArrayPixelSmooth, 1), size(binnedArrayPixelSmooth, 2));
        for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
            for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
                isIN(iPixelX, iPixelY) = inpolygon(projection_view_bins{iProjection}{1}(iPixelX), ...
                    projection_view_bins{iProjection}{2}(iPixelY),...
                    regionLocation(projection_views(iProjection,1),boundary_projection{iProjection}),...
                    regionLocation(projection_views(iProjection,2),boundary_projection{iProjection}));
            end
        end

        figure(1)
        subplot(3, size(regions, 2), iRegion + (size(regions,2) * (iProjection-1)))

        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
        ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
        im = imagesc(projection_view_bins{iProjection}{1}, projection_view_bins{iProjection}{2},...
            binnedArrayPixelSmooth'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        caxis(thisCmap_limits)
        hold on;
    
        %colorbar
        clearvars binnedArrayPixel
        hold on;
        boundary_projection{iProjection} = boundary(regionLocation(projection_views(iProjection,1),:)',...
            regionLocation(projection_views(iProjection,2),:)', 0);
        plot(regionLocation(projection_views(iProjection,1),boundary_projection{iProjection}),...
            regionLocation(projection_views(iProjection,2),boundary_projection{iProjection}),...
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
    
        makepretty;
        clearvars isIN
        caxis(thisCmap_limits)
        set(gca, 'color', [0.5, 0.5, 0.5]);
        xlim([projection_view_bins{iProjection}{1}(1) , projection_view_bins{iProjection}{1}(end) ])
        ylim([projection_view_bins{iProjection}{2}(1) , projection_view_bins{iProjection}{2}(end) ])
        set(gca,'xticklabel',{[]})
        set(gca,'yticklabel',{[]})
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'YColor', [1 1 1])
        set(gca,'XColor', [1 1 1])
   end
  keep passive_data regions thisCmap_limits st av regionResolution structure_alpha theseColors dFR_psth iRegion bregma
end