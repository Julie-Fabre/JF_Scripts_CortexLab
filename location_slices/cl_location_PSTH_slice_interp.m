% get 3 D structure matrx (un smoothed, including NaNs).

function cl_location_PSTH_slice_interp(load_type, figHandle, pcells)

%% cl_locationPSTH
if contains(load_type, 'taskNaive')
    passive_data = load('/home/julie/Dropbox/MATLAB/naive_data5.mat');
    index = 5;
elseif contains(load_type, 'naive')
    passive_data = load('/home/julie/Dropbox/MATLAB/naive_data1.mat');
    index = 1;
elseif contains(load_type, 'taskGo')
    passive_data = load('/home/julie/Dropbox/MATLAB/gogogo_data2.mat');
    index = 2;
elseif contains(load_type, 'taskNoGo')
    passive_data = load('/home/julie/Dropbox/MATLAB/goNogo_data2.mat');
    index = 2;

end

figure(1);
clf;

if nargin < 2 || isempty(figHandle)
    figHandle = figure();
end

if nargin < 3 || isempty(pcells)
    pcells = false;
end

regions = {'CP', 'GPe', 'SNr'};
%pcells = true;
keep passive_data regions pcells index load_type figHandle

figure(figHandle);
clf;

cl_myPaths;
regionResolution = [1, 1, 1, 1, 1, 1, 1];

passive_data.psth_average = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));
%zscore_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
%   ( nanstd(passive_data.psth_average(:, 1:200), [], 2)  +0.001); -> need
%   to zscore all neurons in bin together.


% dFR_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
%     (nanmean(passive_data.psth_average(:, 1:200), 2)  +1);

dFR_psth = (squeeze(passive_data.psth{index}(:, 3, 1, 250:450)) - nanmean(squeeze(passive_data.psth{index}(:, 3, 1, 1:200)), 2)) ./ ...
    (nanmean(squeeze(passive_data.psth{index}(:, 3, 1, 1:200)), 2) + 0.1);

average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, 1, :), 3));
zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);

if pcells
    thisCmap_limits = [-60, 60];
else
    thisCmap_limits = [-20, 20];
end
theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

if ~exist('st', 'var')
    [tv, av, st, bregma] = ya_loadAllenAtlas(atlasBrainRegLocation); 
end

structure_alpha = 0.2;

%% get each region's limits and define chunks
nChunks = 20;
nBins = 30;
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
end

%% AP chunks in ML x DV projection
for iRegion = 1:size(regions, 2)

    curr_plot_structure = st.id(strcmp(st.acronym, regions{iRegion}));

    figure(figHandle);

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

    theseLocations = passive_data.unit_coords;
    theseLocationsBregmaAbs = [(theseLocations(:, 3)), ...
        theseLocations(:, 1), ...
        theseLocations(:, 2)]; % go from AP, DV, ML to ML, AP, DV (like loaded Atlas)

    bregma_ml_point = bregma(3); %2.5 is difference in scaling between
    % brainreg (25 um resolution) and allen (10um resolution, where this bregma value comes from)
    %theseLocationsBregmaAbs=theseLocationsBregmaAbs./2.5;
    theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1) > bregma_ml_point), 1) = ...
        bregma_ml_point - (theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1) > bregma_ml_point), 1))+bregma_ml_point; % squash right hemisphere on the left

    %% plot average increase for each bin

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
                        (passive_data.unitType' == 1); % | passive_data.unitType' ==2);
                else
                    theseNeurons = binX == iBinX & binY == iBinY & ...
                        theseLocationsBregmaAbs(:, 2) >= round(chunks_region(iRegion, iChunk)) & ...
                        theseLocationsBregmaAbs(:, 2) < round(chunks_region(iRegion, iChunk+1)) & ...
                        passive_data.unit_area == iRegion & ...
                        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
                end

                binnedArrayTot = [];
                if sum(theseNeurons) > 0

                    if pcells
                        pcells_2d = sum( ...
                            abs(nanmean(zscore_psth(theseNeurons, 250:350), 2))' > 0.15) / sum(theseNeurons);

                        %      pcells_2d = sum(passive_data.pvalue_shuffled_005{index}(theseNeurons) )/sum(theseNeurons);
                        binnedArrayPixel(iBinX, iBinY) = pcells_2d;
                    else
                        % mean_2d = (nanmean(nanmean(passive_data.psth_average(theseNeurons, 260:360))).*100 - ...
                        %     nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100) ./...
                        %      (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100);%nanmean(abs(dFR_psth(theseNeurons, :)), 2);

                        mean_2d = nanmean(abs((nanmean(passive_data.psth_average(theseNeurons, 260:360), 2) .* 100 - ...
                            nanmean(passive_data.psth_average(theseNeurons, 1:100), 2) .* 100))) ./ ...
                            (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))) .* 100 +0.1);


                        %   mean_2d = (abs(nanmean(passive_data.psth_average(theseNeurons, 260:360))).*100 -  nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100)) ./...
                        %     (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100);%nanmean(abs(dFR_psth(theseNeurons, :)), 2);
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

        binnedArrayPixelSmooth = binnedArrayPixel;

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


        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);

        regionSlices{iRegion}{iChunk} = binnedArrayPixelSmooth;
    end

    % interpolate, smooth and plot
    dataCell = regionSlices{iRegion};
   % nChunks = numel(dataCell);
    dataFilled = dataCell; % Initialize with original data for output

    % Intra-slice Interpolation
    for i = 1:nChunks
        currentSlice = dataCell{i};
        % Using scattered interpolation to handle within-slice NaNs
        [rows, cols] = size(currentSlice);
        [X, Y] = meshgrid(1:cols, 1:rows);
        maskNaN = isnan(currentSlice);
        if any(maskNaN, 'all') % Proceed if there are NaNs
            F = scatteredInterpolant(X(~maskNaN), Y(~maskNaN), currentSlice(~maskNaN), 'linear', 'none');
            currentSlice(maskNaN) = F(X(maskNaN), Y(maskNaN));
        end
        dataFilled{i} = currentSlice; % Store the interpolated slice
    end

    % Attempt to refine NaN filling with information from adjacent slices
    for i = 1:nChunks
        currentSlice = dataFilled{i};
        maskNaN = isnan(currentSlice);
        if any(maskNaN, 'all') % Proceed if there are still NaNs after the initial interpolation
            adjacentValues = [];
            if i > 1 % Look at the previous slice
                prevSlice = dataFilled{i-1};
                adjacentValues = [adjacentValues; prevSlice(maskNaN)]; % Collect corresponding non-NaN values
            end
            if i < nChunks % Look at the next slice
                nextSlice = dataFilled{i+1};
                adjacentValues = [adjacentValues; nextSlice(maskNaN)]; % Collect corresponding non-NaN values
            end

            % Use median or mean of adjacent non-NaN values to fill NaNs
            % Ensuring the dimensions match; you might need to adjust this based on your exact filling strategy
            if ~isempty(adjacentValues)
                medianValues = median(adjacentValues, 2, 'omitnan');
                currentSlice(maskNaN) = medianValues;
                dataFilled{i} = currentSlice; % Update the slice with filled values
            end
        end
    end
    % smooth
    dataSmoothed = cell(nChunks, 1); % Initialize cell array for smoothed data

    % Define Gaussian kernel for smoothing
    % Adjust the size and sigma as per your data's characteristics
    gaussianSize = 5; % Size of the Gaussian filter
    gaussianSigma = 2; % Standard deviation of the Gaussian filter
    G = fspecial('gaussian', gaussianSize, gaussianSigma);

    % Apply smoothing to each slice
    for i = 1:nChunks
        currentSlice = dataFilled{i};
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

    % Assuming each slice has been smoothed within itself as per the previous steps
    % Here, we focus on smoothing across the Z dimension (between slices)

    % Parameters for 1D Gaussian smoothing across slices
    zGaussianSize = 3; % Smaller size as we're only smoothing between a few slices
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

    % Apply 1D Gaussian smoothing across the Z dimension (between slices)
    for x = 1:maxRows
        for y = 1:maxCols
            pixelSeries = squeeze(tempVolume(x, y, :)); % Extract the pixel series across slices
            smoothedPixelSeries = imfilter(pixelSeries, zG, 'same', 'replicate'); % Apply 1D smoothing
            tempVolume(x, y, :) = smoothedPixelSeries; % Update the volume with smoothed values
        end
    end

    % Extract the smoothed slices back into the cell array format
    dataCrossSmoothed = cell(nChunks, 1);
    for i = 1:nChunks
        dataCrossSmoothed{i} = tempVolume(:, :, i);
    end


    %     figure();
    % for iChunk=1:nChunks
    %     subplot(1,nChunks,iChunk)
    %      imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
    %             regionSlices{iRegion}{iChunk}'.*100);
    %        % set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    %
    %        % set(gca, 'color', [0.5, 0.5, 0.5]);
    %         colormap(brewermap([], '*RdBu'));
    %         caxis(thisCmap_limits)
    % end
    %
    % figure();
    % for iChunk=1:nChunks
    %     subplot(1,nChunks,iChunk)
    %      imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
    %             dataFilled{iChunk}'.*100);
    %        % set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    %
    %        % set(gca, 'color', [0.5, 0.5, 0.5]);
    %         colormap(brewermap([], '*RdBu'));
    %         caxis(thisCmap_limits)
    % end
    %
    % figure();
    % for iChunk=1:nChunks
    %     subplot(1,nChunks,iChunk)
    %      imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
    %             dataSmoothed{iChunk}'.*100);
    %        % set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    %
    %        % set(gca, 'color', [0.5, 0.5, 0.5]);
    %         colormap(brewermap([], '*RdBu'));
    %         caxis(thisCmap_limits)
    % end
    %
    %    figure();
    % for iChunk=1:nChunks
    %     subplot(1,nChunks,iChunk)
    %      imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
    %             dataCrossSmoothed{iChunk}'.*100);
    %        % set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    %
    %        % set(gca, 'color', [0.5, 0.5, 0.5]);
    %         colormap(brewermap([], '*RdBu'));
    %         caxis(thisCmap_limits)
    % end

    dataCrossSmoothed_regions{iRegion} = dataCrossSmoothed;

    %for iRegion =1:3
    for iChunk = 1:nChunks
        %else
        %    binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
        %end
        subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
        ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            dataCrossSmoothed_regions{iRegion}{iChunk}'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));

        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        caxis(thisCmap_limits)
        hold on;

        %colorbar
        clearvars binnedArrayPixel
        hold on;
         region_area = permute(av(round(chunks_region(iRegion, iChunk)):1:round(chunks_region(iRegion, iChunk+1)), ...
            1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
        % AP, DV, ML -> ML, AP, DV
clearvars regionLocation
        [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
            = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
 

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

        axis off;

        % Alternatively, to remove only specific components:
        set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
        xlabel(''); % Remove x-axis label
        ylabel(''); % Remove y-axis label
        % title(''); % Remove title

        %end

    end
     keep dataCrossSmoothed_regions figHandle passive_data regions thisCmap_limits st av regionResolution structure_alpha theseColors dFR_psth iRegion bregma nChunks chunks_region pcells index load_type zscore_psth
  
    prettify_plot
end