
%% cl_locationPSTH
if contains(load_type, 'naive')
    passive_data = load('/home/julie/Dropbox/MATLAB/naive_data1.mat');
    index = 1;
elseif contains(load_type, 'taskGo')
    passive_data = load('/home/julie/Dropbox/MATLAB/gogogo_data2.mat');
    index = 2;
elseif contains(load_type, 'taskNoGo')
    passive_data = load('/home/julie/Dropbox/MATLAB/goNogo_data2.mat');
    index = 2;
end

pcells = true;
keep passive_data pcells index load_type
figure(1);
clf;
figure(2);
clf;

cl_myPaths;
regionResolution = [1, 1, 1, 1, 1, 1, 1];

passive_data.psth_average = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));
%zscore_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
%   ( nanstd(passive_data.psth_average(:, 1:200), [], 2)  +0.001); -> need
%   to zscore all neurons in bin together.
dFR_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
    (nanmean(passive_data.psth_average(:, 1:200), 2) + 0.001);
if pcells
    thisCmap_limits = [-100, 100];
else
    thisCmap_limits = [-65, 65];
end
theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

if ~exist('st', 'var')
    [tv, av, st, bregma] = ya_loadAllenAtlas(atlasBrainRegLocation); 
end

structure_alpha = 0.2;


cl_myPaths;

%% get loading info
% recording type exp info
if contains(load_type, 'naive')
    info_table = readtable([csvPath, 'allPassiveRecs.csv'], 'VariableNamingRule', 'modify');
elseif contains(load_type, 'taskGo')

    info_table = readtable([csvPath, 'allTaskGoGoGoRecs.csv'], 'VariableNamingRule', 'modify');
elseif contains(load_type, 'taskNoGo')
    info_table = readtable([csvPath, 'allTaskRecs.csv'], 'VariableNamingRule', 'modify');
end

% regions to load
regions = {'CP', 'GPe', 'SNr'}; % 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
warning off;
regions_id = [672, 1022, 381];

% get all mice x thisDate x site combinations
use_recs_full = strcmp(info_table.Use_, 'Yes');

info_table = sortrows(info_table, 1); % make sure data is sorted by mouse
unique_mice = unique(info_table.Mouse(use_recs_full), 'stable');

%% get each region's limits and define chunks
nChunks = 10;
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
region_ap_boundaries = [152, 307; ...
    219, 287; ...
    311, 363];

for iRegion = 1:size(regions, 2)
    chunks_region(iRegion, :) = region_ap_boundaries(iRegion, 1):(region_ap_boundaries(iRegion, 2) - region_ap_boundaries(iRegion, 1)) / nChunks:region_ap_boundaries(iRegion, 2);
end

%% Visualize AP chunks in ML x AP projection

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

        subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (nChunks - iChunk)))
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

    prettify_plot('XLimits', 'col');
end

%% AP chunks in ML x DV projection

for iMouse = 3%1:7%:size(unique_mice,1)
    mouse = unique_mice{iMouse};
    curr_fig = figure('Name', mouse);
    for iRegion = 1:size(regions, 2)

        curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));


        projection_views = repmat([1, 3], nChunks, 1); %[1,2; 1,3; 2,3];%ML/AP, ML/DV, AP/DV

        figure(curr_fig)
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

            subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (nChunks - iChunk)))
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
        %theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
        %    theseLocations(:, 1), ...
        %    theseLocations(:, 2)]; % go from AP, DV, ML to ML, AP, DV (like loaded Atlas)

        theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
            theseLocations(:, 1), ...
            theseLocations(:, 2)]; % go from AP, DV, ML to ML, AP, DV (like loaded Atlas)
        
        bregma_ml_point = bregma(3); %2.5 is difference in scaling between
        % brainreg (25 um resolution) and allen (10um resolution, where this bregma value comes from)

        theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1)>bregma_ml_point),1) = ...
            bregma_ml_point - ...
            (theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1)>bregma_ml_point), 1))+...
            + bregma_ml_point;% squash right hemisphere on the left



        
        %% plot average increase for each bin
       % try
            for iChunk = 1:nChunks
                [N, Xedges, Yedges, binX, binY] = histcounts2(theseLocationsBregmaAbs(:, projection_views(iChunk, 1)), ...
                    theseLocationsBregmaAbs(:, projection_views(iChunk, 2)), projection_view_bins{iChunk}{1}, ...
                    projection_view_bins{iChunk}{2}); %par rapport a bregma!


                binnedArrayPixel = nan(size(Xedges, 2), size(Yedges, 2)); % initialize

                for iBinX = 1:size(Xedges, 2)
                    for iBinY = 1:size(Yedges, 2)
                        if pcells
                            theseNeurons = binX == iBinX & binY == iBinY & ...
                                theseLocationsBregmaAbs(:, 2) >= round(chunks_region(iRegion, iChunk)) & ...
                                theseLocationsBregmaAbs(:, 2) < round(chunks_region(iRegion, iChunk+1)) & ...
                                passive_data.unit_area == iRegion & ...
                                (passive_data.unitType' == 1) & passive_data.animal_thisDate_site_shank(:, 1) == iMouse; % | passive_data.unitType' ==2);
                        else
                            theseNeurons = binX == iBinX & binY == iBinY & ...
                                theseLocationsBregmaAbs(:, 2) >= round(chunks_region(iRegion, iChunk)) & ...
                                theseLocationsBregmaAbs(:, 2) < round(chunks_region(iRegion, iChunk+1)) & ...
                                passive_data.unit_area == iRegion & ...
                                (passive_data.unitType' == 1 | passive_data.unitType' == 2) & passive_data.animal_thisDate_site_shank(:, 1) == iMouse;
                        end

                        binnedArrayTot = [];
                        if sum(theseNeurons) > 0

                            if pcells
                                pcells_2d = sum(passive_data.pvalue_shuffled_005{index}(theseNeurons)) / sum(theseNeurons);
                                binnedArrayPixel(iBinX, iBinY) = pcells_2d;
                            else
                                mean_2d = nanmean(abs(dFR_psth(theseNeurons, :)), 2);
                                binnedArrayPixel(iBinX, iBinY) = nanmean(mean_2d(~isinf(mean_2d)));
                            end

                        end

                    end
                end

                % %NaN-out 0 values (= no data) and Inf values if present (QQ need to check why Inf values sometimes crop up)
                % binnedArrayPixel(binnedArrayPixel == Inf) = NaN;
                if ~pcells
                    binnedArrayPixel(binnedArrayPixel == 0) = NaN;
                end

                % smooth data
                binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 2, 2);
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

                figure(curr_fig);
                subplot(nChunks, size(regions, 2), iRegion+(size(regions, 2) * (nChunks - iChunk)))
                if pcells
                    binnedArrayPixelSmooth(isIN == 0) = min(thisCmap_limits);
                else
                    binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
                end
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
        %catch
        %end
    end
    keep passive_data regions thisCmap_limits st av regionResolution structure_alpha theseColors dFR_psth iRegion bregma nChunks chunks_region pcells index unique_mice load_type
end
  