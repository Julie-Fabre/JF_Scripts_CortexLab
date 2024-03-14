% get 3 D structure matrx (un smoothed, including NaNs).

function cl_location_PSTH_slice_interp_nans(load_type, figHandle, pcells, incDec)
%QQ add - use no move trials!

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
elseif strcmp(load_type, 'task')
    passive_data = load('/home/julie/Dropbox/MATLAB/gogogo_data2.mat');
    passive_data2 = load('/home/julie/Dropbox/MATLAB/goNogo_data2.mat');
    index = 2;
    passive_data.unitType = [passive_data.unitType, passive_data2.unitType];
    passive_data.unit_area = [passive_data.unit_area; passive_data2.unit_area];
    passive_data.unit_coords = [passive_data.unit_coords; passive_data2.unit_coords];
    passive_data.psth{index} = [passive_data.psth{index}; passive_data2.psth{index}];
    passive_data.pvalue_shuffled_005{index} = [passive_data.pvalue_shuffled_005{index}, passive_data2.pvalue_shuffled_005{index}];
end

figure(1);
clf;

if nargin < 2 || isempty(figHandle)
    figHandle = figure();
end

if nargin < 3 || isempty(pcells)
    pcells = false;
end

if nargin < 4 || isempty(incDec)
    incDec = false;
end

regions = {'CP', 'GPe', 'SNr'};
%pcells = true;
keep passive_data regions pcells index load_type figHandle incDec

figure(figHandle);
figHandle.Color = 'w';
clf;

cl_myPaths;
regionResolution = [1, 1, 1, 1, 1, 1, 1];

passive_data.psth_average = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));
%zscore_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
%   ( nanstd(passive_data.psth_average(:, 1:200), [], 2)  +0.001); -> need
%   to zscore all neurons in bin together.


% dFR_psth = (passive_data.psth_average(:, 250:450) - nanmean(passive_data.psth_average(:, 1:200), 2)) ./ ...
%     (nanmean(passive_data.psth_average(:, 1:200), 2)  +1);

%dFR_psth = (squeeze(passive_data.psth{index}(:, 3, 1, 250:450)) - nanmean(squeeze(passive_data.psth{index}(:, 3, 1, 1:200)), 2)) ./ ...
%    (nanmean(squeeze(passive_data.psth{index}(:, 3, 1, 1:200)), 2) + 0.1);

average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, 1, :), 3));
zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);

if pcells
    thisCmap_limits = [-30, 30];
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
    chunks_region_ap(iRegion, :) = (bregma(1) - chunks_region(iRegion, :)) * 25;
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

for iRegion = 1:size(regions, 2)

    %% get data
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
        bregma_ml_point - (theseLocationsBregmaAbs(find(theseLocationsBregmaAbs(:, 1) > bregma_ml_point), 1)) + bregma_ml_point; % squash right hemisphere on the left

    %% interpolate and smooth
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
                        if sum(theseNeurons) > 1
                            pcells_2d = sum( ...
                                abs(nanmean(zscore_psth(theseNeurons, 250:350), 2))' > 0.25 & ...
                                passive_data.pvalue_shuffled_005{index}(theseNeurons)) / sum(theseNeurons);
                        else
                            pcells_2d = 0;
                        end

                        %      pcells_2d = sum(passive_data.pvalue_shuffled_005{index}(theseNeurons) )/sum(theseNeurons);
                        binnedArrayPixel(iBinX, iBinY) = pcells_2d;
                    elseif incDec
                        mean_2d_inc = nanmean((nanmean(passive_data.psth_average(theseNeurons, 260:360), 2) .* 100 - ...
                            nanmean(passive_data.psth_average(theseNeurons, 1:100), 2) .* 100)) ./ ...
                            (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))) .* 100 + 0.1);
                    else
                        % mean_2d = (nanmean(nanmean(passive_data.psth_average(theseNeurons, 260:360))).*100 - ...
                        %     nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100) ./...
                        %      (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))).*100);%nanmean(abs(dFR_psth(theseNeurons, :)), 2);

                        mean_2d = nanmean(abs((nanmean(passive_data.psth_average(theseNeurons, 260:360), 2) .* 100 - ...
                            nanmean(passive_data.psth_average(theseNeurons, 1:100), 2) .* 100))) ./ ...
                            (nanmean(nanmean(passive_data.psth_average(theseNeurons, 1:100))) .* 100 + 0.1);


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


        binnedArrayPixelSmooth(isIN == 0) = Inf;

        regionSlices{iRegion}{iChunk} = binnedArrayPixelSmooth;
    end

    % interpolate, smooth and plot
    dataCell = regionSlices{iRegion};
    % nChunks = numel(dataCell);
    dataFilled = dataCell; % Initialize with original data for output

    if iRegion==2 &&  contains(load_type, 'naive')%temporary really crappy Qq remove!!!
        dataCell{15}(find(~isnan(dataFilled{15}) & ~isinf(dataFilled{15})))=NaN;
         dataCell{16}(find(~isnan(dataFilled{16}) & ~isinf(dataFilled{16})))=NaN;
        dataFilled{15}(find(~isnan(dataFilled{15}) & ~isinf(dataFilled{15})))=NaN;
        dataFilled{16}(find(~isnan(dataFilled{15}) & ~isinf(dataFilled{15})))=NaN;
    end

    % Intra-slice Interpolation
    for i = 1:nChunks
        currentSlice = dataCell{i};
        % Using scattered interpolation to handle within-slice NaNs
        [rows, cols] = size(currentSlice);
        [X, Y] = meshgrid(1:cols, 1:rows);
        maskNaN = ~isnan(currentSlice) & ~isinf(currentSlice);
        stayIntheseLims = isinf(currentSlice);
        maskNaN_exists = isnan(currentSlice);
        if any(maskNaN_exists, 'all') && any(maskNaN, 'all') % if not only NaNs/Infs, otherwise, skip
            se = strel('square', 3); % A 3x3 square will dilate the mask by 1 pixel all around
            maskNaN_dilated = imdilate(maskNaN, se);
            maskNaN_dilated(stayIntheseLims) = 0;
            try
                F = scatteredInterpolant(X(maskNaN), Y(maskNaN), currentSlice(maskNaN), 'linear', 'none');
                currentSlice(maskNaN_dilated & ~maskNaN) = F(X(maskNaN_dilated & ~maskNaN), Y(maskNaN_dilated & ~maskNaN));
            catch
            end
            dataFilled{i} = currentSlice; % Store the interpolated slice

        end
        %maskNaN = isnan(currentSlice)
        % if any(maskNaN, 'all') % Proceed if there are NaNs
        %     F = scatteredInterpolant(X(~maskNaN), Y(~maskNaN), currentSlice(~maskNaN), 'linear', 'none');
        %     currentSlice(maskNaN) = F(X(maskNaN), Y(maskNaN));
        % end
        dataFilled{i} = currentSlice; % Store the interpolated slice
    end

    % Attempt to refine NaN filling with information from adjacent slices
    for i = 1:nChunks
        currentSlice = dataFilled{i};
        %maskNaN = isnan(currentSlice);
        maskNaN_exists = isnan(currentSlice);
        maskNaN = ~isnan(currentSlice) & ~isinf(currentSlice);
        if any(maskNaN_exists, 'all') && any(maskNaN, 'all') % if not only NaNs/Infs, otherwise, skip
            adjacentValues = [];
            se = strel('square', 9); % A 3x3 square will dilate the mask by 1 pixel all around
            maskNaN_dilated = imdilate(maskNaN, se);
            maskNaN_dilated(stayIntheseLims) = 0;
            if i > 1 % Look at the previous slice
                prevSlice = dataFilled{i-1};
                adjacentValues = [adjacentValues, prevSlice(maskNaN_dilated & ~maskNaN)]; % Collect corresponding non-NaN values
            end
            if i < nChunks % Look at the next slice
                nextSlice = dataFilled{i+1};
                adjacentValues = [adjacentValues, nextSlice(maskNaN_dilated & ~maskNaN)]; % Collect corresponding non-NaN values
            end

            % Use median or mean of adjacent non-NaN values to fill NaNs
            % Ensuring the dimensions match; you might need to adjust this based on your exact filling strategy

            if ~isempty(adjacentValues)
                maskNaN_donut = maskNaN_dilated & ~maskNaN;
                adjacentValues(isinf(adjacentValues)) = NaN;
                medianValues = median(adjacentValues, 2, 'omitnan');
                currentSlice(maskNaN_donut) = medianValues;
                dataFilled{i} = currentSlice; % Update the slice with filled values
            end
        end
    end
    % smooth
    dataSmoothed = cell(nChunks, 1); % Initialize cell array for smoothed data

    % Define Gaussian kernel for smoothing
    % Adjust the size and sigma as per your data's characteristics
    if iRegion==1
    gaussianSize = 5; % Size of the Gaussian filter
    gaussianSigma = 2; % Standard deviation of the Gaussian filter
    else
        gaussianSize = 5; % Size of the Gaussian filter
        gaussianSigma = 4; % Standard deviation of the Gaussian filter
    end
    G = fspecial('gaussian', gaussianSize, gaussianSigma);

    % Apply smoothing to each slice
    for i = 1:nChunks
        currentSlice = dataFilled{i};
        if ~isempty(currentSlice) % Ensure the slice is not empty
            % Apply Gaussian smoothing
            currentSlice(isinf(currentSlice)) = NaN;
            %smoothedSlice = imfilter(currentSlice, G, 'same', 'replicate');
            smoothedSlice = nanconv(currentSlice, G, 'edge', 'nanout');

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
    if iRegion==1
    zGaussianSize = 3; % Smaller size as we're only smoothing between a few slices
    zGaussianSigma = 1; % Adjust based on the desired level of cross-slice smoothing
    else
         zGaussianSize = 4; % Smaller size as we're only smoothing between a few slices
    zGaussianSigma = 2; % Adjust based on the desired level of cross-slice smoothing
    end
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
            %smoothedPixelSeries = imfilter(pixelSeries, zG, 'same', 'replicate'); % Apply 1D smoothing
            smoothedPixelSeries = smooth1DWithNaN(pixelSeries, zG); % Apply custom 1D smoothing with NaN handling
            tempVolume(x, y, :) = smoothedPixelSeries; % Update the volume with smoothed values
        end
    end

    % Extract the smoothed slices back into the cell array format
    dataCrossSmoothed = cell(nChunks, 1);
    for i = 1:nChunks
        dataCrossSmoothed{i} = tempVolume(:, :, i);
    end

    dataCrossSmoothed_regions{iRegion} = dataCrossSmoothed;

    %% plot data
    if iRegion > 1 && ~pcells
        thisCmap_limits = [-12, 12];
    end
    for iChunk = 1:nChunks
        clearvars regionLocation
        % get structure boundaries and plot outline
        region_area = permute(av(round(chunks_region(iRegion, iChunk)):1:round(chunks_region(iRegion, iChunk+1)), ...
            1:1:end, 1:1:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere


        [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
            = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
        isIN = nan(size(dataCrossSmoothed_regions{iRegion}{iChunk}, 1), size(dataCrossSmoothed_regions{iRegion}{iChunk}, 2));
        for iPixelX = 1:size(dataCrossSmoothed_regions{iRegion}{iChunk}, 1)
            for iPixelY = 1:size(dataCrossSmoothed_regions{iRegion}{iChunk}, 2)
                isIN(iPixelX, iPixelY) = inpolygon(projection_view_bins{iChunk}{1}(iPixelX), ...
                    projection_view_bins{iChunk}{2}(iPixelY), ...
                    regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
                    regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}));
            end
        end


        dataCrossSmoothed_regions{iRegion}{iChunk}(isIN == 0) = 0;

        %else
        %    binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap_limits);
        %end
        subplot(size(regions, 2), nChunks, (iRegion - 1).*nChunks+iChunk)
        ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red

        % Create a gray background layer
        bg = ones(size(dataCrossSmoothed_regions{iRegion}{iChunk})) .* -Inf; % 0.5 corresponds to mid-gray in many colormaps
        imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, bg); % Plot the gray background
        hold on; % Keep the gray background, and plot on top of it

        %projection_view_bins{iChunk}(isinf(projection_view_bins{iChunk})) =0;
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            dataCrossSmoothed_regions{iRegion}{iChunk}'*100);

        % NaN values transparent
        %set(gcf, 'color', [0.5, 0.5, 0.5]); -> this sets whole figure to
        %gray.
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));

        % Generate the colormap using brewermap
        originalColormap = brewermap([], '*RdBu');
        middleIndex = ceil(size(originalColormap, 1)/2);
        whiteColor = [1, 1, 1];
        grayColor = [0.9, 0.9, 0.9];
        modifiedColormap = [grayColor; originalColormap(1:middleIndex-1, :); whiteColor; whiteColor; originalColormap(middleIndex+1:end, :)];

        % Apply the modified colormap
        colormap(modifiedColormap);

        % Adjust caxis if necessary to include the extra colors added
        % Assume thisCmap_limits is defined for your data's range
        caxis([min(thisCmap_limits), max(thisCmap_limits)]);


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
            'Color', theseColors{iRegion}, 'LineWidth',2);

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
        %caxis(thisCmap_limits)
        %set(gca, 'color', [0.5, 0.5, 0.5]);
        xlim([projection_view_bins{iChunk}{1}(1), projection_view_bins{iChunk}{1}(end)])
        ylim([projection_view_bins{iChunk}{2}(1), projection_view_bins{iChunk}{2}(end)])
        set(gca, 'Layer', 'top');
        axis off;

        % Alternatively, to remove only specific components:
        set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
        xlabel(''); % Remove x-axis label
        ylabel(''); % Remove y-axis label
        thisAp_slice = round(nanmean(chunks_region_ap(iRegion, iChunk:iChunk+1))./1000, 2);
        if iChunk == 1
            title(['Bregma: ', num2str(thisAp_slice)]); % Remove title
        elseif iChunk == nChunks
            title([num2str(thisAp_slice), 'mm']); % Remove title
        else
            title([num2str(thisAp_slice)]); % Remove title
        end


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

    keep incDec chunks_region_ap dataCrossSmoothed_regions figHandle passive_data regions thisCmap_limits st av regionResolution structure_alpha theseColors iRegion bregma nChunks chunks_region pcells index load_type zscore_psth

    %prettify_plot
end


end