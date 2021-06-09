keep ephysData ephysDataSnr

animalsType = {'Naive'};
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
regionSpacing = [20, 5, 5, 5, 5];
%regions = {'DMS', 'PS', 'STN', 'GPe', 'SNr', 'GPi'};
recordingInfo = readtable('C:\Users\Julie\Dropbox\Analysis\Recordings - Sheet1.csv');
bregma = [540, 0, 570];
thisCmap = [-20, 20];

%% X-Y, X-Z, Y-Z plots for each region. bin by 100 um (?) and plot in colormap increase in response ?

for iType = 1%:size(animalsType, 2)
    theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

    allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
    tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
    av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
    slice_spacing = 10;
    structure_alpha = 0.2;
    %figure();
    %overlay regions
    smoothF = 3;
    for iRegion = 1:size(regions, 2)
        if iRegion ==1
            thisData = ephysData;%.ephysData;
        elseif iRegion ==4
            thisData = ephysDataSnr;
        else
            thisData=ephysData;
        end
        curr_plot_structure = find(strcmp(st.acronym, regions{iRegion}));
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
        %         plotBrainGrid();
        %         hold on;
        %         scatter3(bregma(1), bregma(3), bregma(2), 'filled')
        %         hold on;
        %         axis vis3d equal off manual
        %         view([-30, 25]);
        %         caxis([0, 300]);
        %         [ap_max, dv_max, ml_max] = size(tv);
        %         xlim([-10, ap_max + 10])
        %         ylim([-10, ml_max + 10])
        %         zlim([-10, dv_max + 10])
        %         structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
        %             'Faces', structure_3d.faces, ...
        %             'FaceColor', theseColors{iRegion, :}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
        %

        % 2-D
        ii = permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
        [r, c, v] = ind2sub(size(ii), find(ii));

        figure(1);
        % x-y
        subplot(3, size(regions, 2), iRegion)
        bb = boundary(r, c, 0);
        plot(r(bb)./10-(bregma(3) / 100), c(bb)./10-(bregma(1) / 100), 'Color', theseColors{iRegion});
        %xlabel('ML (mm)')
        %ylabel('AP (mm)')
        axis equal
        axis square
        axis image


        xyXLim = xlim;
        xyYLim = ylim;
        makepretty;

        % x-z
        subplot(3, size(regions, 2), size(regions, 2)+iRegion)
        bb2 = boundary(r, v, 0);
        plot(r(bb2)./10-(bregma(3) / 100), -v(bb2)./10, 'Color', theseColors{iRegion});
        %xlabel('ML (mm)')
        %ylabel('DV (mm)')
        axis equal
        axis square
        axis image
        makepretty;
        xzXLim = xlim;
        xzZLim = ylim;

        % y-z
        subplot(3, size(regions, 2), size(regions, 2)*2+(iRegion))
        bb3 = boundary(c, v, 0);
        plot(c(bb3)./10-(bregma(1) / 100), -v(bb3)./10, 'Color', theseColors{iRegion});
        %xlabel('AP (mm)')
        %ylabel('DV (mm)')
        axis equal
        axis square
        axis image
        makepretty;
        yzYLim = xlim;
        yzZLim = ylim;
        % bin region in 200um bins
        xy = [r(bb), c(bb)];

        xyCountBins = {xyXLim(1) * 10 :(xyXLim(2) - xyXLim(1)) * 10 / 15:xyXLim(2) * 10, xyYLim(1) * 10 :(xyYLim(2) - xyYLim(1)) * 10 / 15:xyYLim(2) * 10};
        xzCountBins = {xzXLim(1) * 10 :(xzXLim(2) - xzXLim(1)) * 10 / 15:xzXLim(2) * 10, xzZLim(1) * 10 :(xzZLim(2) - xzZLim(1)) * 10 / 15:xzZLim(2) * 10};
        yzCountBins = {yzYLim(1) * 10 :(yzYLim(2) - yzYLim(1)) * 10 / 15:yzYLim(2) * 10, yzZLim(1) * 10 :(yzZLim(2) - yzZLim(1)) * 10 / 15:yzZLim(2) * 10};
        % load all region data [-0.1 - 0.3]s after stim onset, with which bin belongs to

        CNew = cat(1, thisData.location);

        ind = find(ismember(CNew, regions{iRegion}));
       % clearvars theseLocationsInfo
        if iRegion == 4
        for iLocation = 1:size(ind, 1)
            theseLocationsInfo(iLocation) = size(thisData(ind(iLocation)).template_location, 1)+1;
        end
        else
        for iLocation = 1:size(ind, 1)
            theseLocationsInfo(iLocation) = size(thisData(ind(iLocation)).template_location, 1);
        end
        end
        theseLocationsInfo = cumsum(theseLocationsInfo);
        theseLocations = cat(1, thisData(ind).template_location);
        theseLocationsBregmaAbs = [abs(theseLocations(:, 3)-bregma(3)), abs(theseLocations(:, 1)-bregma(1)), theseLocations(:, 2)];


        [N, Xedges, Yedges, binX, binY] = histcounts2(((bregma(3))-abs(theseLocations(:, 3)-(bregma(3)))-bregma(3))./100, (theseLocations(:, 1)-bregma(1))./100, xyCountBins{1, 1}./10, xyCountBins{1, 2}./10); %par rapport a bregma!
        thisSpacer = 0;
        for iBinX = 1:size(Xedges, 2)
            for iBinY = 1:size(Yedges, 2)
                theseNeurons = binX == iBinX & binY == iBinY;
                theseNeuronsInd = find(theseNeurons);
                binnedArrayTot = [];
                if ~isempty(theseNeuronsInd)
                    [h, hh] = histc(theseNeuronsInd, theseLocationsInfo);
                    uniqueRecs = unique(hh) + 1;
                    for iUniqueRec = 1:size(uniqueRecs, 1) %get psth per rec
                        theseTheseNeurons = theseNeuronsInd(hh == uniqueRecs(iUniqueRec));
<<<<<<< HEAD
                        if ~isempty(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(thisData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(thisData(uniqueRecs(iUniqueRec)).spike_templates, ...
=======
                        if ~isempty(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(ephysData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(ephysData(uniqueRecs(iUniqueRec)).spike_templates, ...
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec))-1);
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                thisData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                            binnedArrayTot = [binnedArrayTot; binnedArray];
                            
                            hold on;
                        end
                    end

                end
                if ~isempty(binnedArrayTot)
                    binnedArrayTot = binnedArrayTot(any(binnedArrayTot,2),:);
                    binnedArrayTotBinned(iBinX, iBinY, :) = nanmean(binnedArrayTot, 1); %average
                    figure(1+iRegion);
                    subplot(311)
                    if isnan(thisSpacer)
                        thisSpacer = 0;
                    end
                    plot(-0.2:0.01:0.3-0.01, nanmean(binnedArrayTot, 1)+thisSpacer )
                    thisSpacer = thisSpacer + nanmax(nanmean(binnedArrayTot, 1));
                            
                    % figure();
                    % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinY, :)))
                    % if all(nanmean(binnedArrayTotBinned(iBinX, iBinY, :)))
                    % figure();
                    % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinZ, :)))
                    binnedArrayPixel(iBinX, iBinY) = (nanmean(binnedArrayTotBinned(iBinX, iBinY, 26:41)) - ...
                        nanmean(binnedArrayTotBinned(iBinX, iBinY, 1:20))) ./ nanmean(binnedArrayTotBinned(iBinX, iBinY, 1:20));
                    
                    %else
                    %     binnedArrayTotBinned(iBinX, iBinY, :) = nan(50, 1); %average
                    %binnedArrayPixel(iBinX, iBinY) = NaN;
                    %end
                else
                    binnedArrayTotBinned(iBinX, iBinY, :) = nan(50, 1); %average
                    binnedArrayPixel(iBinX, iBinY) = NaN;
                end
            end

        end
%         figure(iRegion+1)
%         subplot(311)
%         %plot(squeeze(binnedArrayTotBinned(iBinX, iBinY,:)))
%         rrr = reshape(binnedArrayTotBinned(:, :,:), size(binnedArrayTotBinned,1).*size(binnedArrayTotBinned,2),size(binnedArrayTotBinned,3));
%         rrrr = (rrr(any(rrr,2),:)-nanmean(rrr(any(rrr,2),26:41),2))./nanmean(rrr(any(rrr,2),26:41),2);
%         imagesc(-0.2:0.01:0.3-0.01,[],rrrr)    
%         colormap(brewermap([], '*RdBu'));
%         caxis([-max(max(abs(rrrr))) max(max(abs(rrrr)))])
        %subplot(212)
        %plot(-0.2:0.01:0.3-0.01, nanmean(squeeze(nanmean(binnedArrayTotBinned(:, :, :)))), 'k', 'LineWidth', 2)
        figure(1)

        binnedArrayPixel(binnedArrayPixel == Inf) = NaN;
<<<<<<< HEAD
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, smoothF, smoothF);
=======
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 4, 4);
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3

        subplot(3, size(regions, 2), iRegion)

        rbb = r(bb);
        cbb = c(bb);
        for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
            for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
                %                 thisPixelLocation = [xyCountBins{1, 1}(iPixelX)/10, xyCountBins{1, 2}(iPixelY)/10];
                %                 [minDistance, indexOfMin] = min(abs(rbb-thisPixelLocation));
                %                % for itt = 1:length(minDistance
                %                indexOfMin = unique(indexOfMin);
                %                thisV = rbb(indexOfMin);
                %                theseVfind(rbb==rbb(indexOfMin);
                %                for iV = 1:length(theseVfind)
                %                end
                %                [~, v2] = min(abs(rbb([1:indexOfMin-1,indexOfMin+1:end])-thisV));
                %                 XatThisY = cbb(rbb(indexOfMin));
                %
                %                  isINX = thisPixelLocation(2)
                %
                isIN(iPixelX, iPixelY) = inpolygon(xyCountBins{1, 1}(iPixelX)./10, ...
                    xyCountBins{1, 2}(iPixelY)./10, r(bb)./10-(bregma(3) / 100), c(bb)./10-(bregma(1) / 100));
            end
        end

        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap);
                ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
<<<<<<< HEAD
        box off
      %  set(
=======
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        im = imagesc(xyCountBins{1, 1}./10, xyCountBins{1, 2}./10, binnedArrayPixelSmooth'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));

        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        caxis(thisCmap)
        hold on;

        %colorbar
        clearvars binnedArrayPixel
        hold on;
        bb = boundary(r, c, 0);
        plot(r(bb)./10-(bregma(3) / 100), c(bb)./10-(bregma(1) / 100), 'Color', theseColors{iRegion});
        xlabel('ML (mm)')
        ylabel('AP (mm)')
        axis equal
        axis square
        axis image

<<<<<<< HEAD
%         nColors = numel(ax.YTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
%         end

%         nColors = numel(ax.XTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
%         end
%         ax.XLabel.Color = [0, 0, 0];
%         ax.YLabel.Color = [0, 0, 0];
=======
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
        ax.XLabel.Color = [0, 0, 0];
        ax.YLabel.Color = [0, 0, 0];
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        makepretty;
        clearvars isIN
        caxis(thisCmap)
        set(gca, 'color', [0.5, 0.5, 0.5]);
        xlim([xyCountBins{1, 1}(1) ./ 10, xyCountBins{1, 1}(end) ./ 10])
        ylim([xyCountBins{1, 2}(1) ./ 10, xyCountBins{1, 2}(end) ./ 10])


        [N, Xedges, Zedges, binX, binZ] = histcounts2(((bregma(3))-abs(theseLocations(:, 3)-(bregma(3)))-bregma(3))./100, -theseLocations(:, 2)./100, xzCountBins{1, 1}./10, xzCountBins{1, 2}./10); %par rapport a bregma!

        for iBinX = 1:size(Xedges, 2)
            for iBinZ = 1:size(Zedges, 2)
                theseNeurons = binX == iBinX & binZ == iBinZ;
                theseNeuronsInd = find(theseNeurons);
                binnedArrayTot = [];
                if ~isempty(theseNeuronsInd)
                    [h, hh] = histc(theseNeuronsInd, theseLocationsInfo);
                    uniqueRecs = unique(hh) + 1;
                    for iUniqueRec = 1:size(uniqueRecs, 1) %get psth per rec
                        theseTheseNeurons = theseNeuronsInd(hh == uniqueRecs(iUniqueRec));
                        if ~isempty(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(thisData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(thisData(uniqueRecs(iUniqueRec)).spike_templates, ...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec))-1);
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                thisData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                            binnedArrayTot = [binnedArrayTot; binnedArray];
                            %figure();
                            %                    plot(-0.2:0.01:0.3-0.01,psth)
                        end
                    end

                end
                if ~isempty(binnedArrayTot)
                    binnedArrayTotBinned(iBinX, iBinZ, :) = nanmean(binnedArrayTot, 1); %average
                    if all(nanmean(binnedArrayTotBinned(iBinX, iBinZ, :)))
                        % figure();
                        % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinZ, :)))
                        binnedArrayPixel(iBinX, iBinZ) = (nanmean(binnedArrayTotBinned(iBinX, iBinZ, 26:41)) - ...
                            nanmean(binnedArrayTotBinned(iBinX, iBinZ, 1:20))) ./ nanmean(binnedArrayTotBinned(iBinX, iBinZ, 1:20));
                    else
                        binnedArrayTotBinned(iBinX, iBinZ, :) = nan(50, 1); %average
                        binnedArrayPixel(iBinX, iBinZ) = NaN;
                    end

                else
                    binnedArrayTotBinned(iBinX, iBinZ, :) = nan(50, 1); %average
                    binnedArrayPixel(iBinX, iBinZ) = NaN;
                end
            end
        end
<<<<<<< HEAD
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, smoothF, smoothF);
=======
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 4, 4);
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        subplot(3, size(regions, 2), size(regions, 2)+iRegion)
        for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
            for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
                isIN(iPixelX, iPixelY) = inpolygon(xzCountBins{1, 1}(iPixelX), ...
                    xzCountBins{1, 2}(iPixelY), r(bb2)-(bregma(3) / 10), -v(bb2));
            end
        end
         figure(iRegion+1)
        subplot(312)
        %plot(squeeze(binnedArrayTotBinned(iBinX, iBinY,:)))
        rrr = reshape(binnedArrayTotBinned(:, :,:), size(binnedArrayTotBinned,1).*size(binnedArrayTotBinned,2),size(binnedArrayTotBinned,3));
        rrrr = (rrr(any(rrr,2),:)-nanmean(rrr(any(rrr,2),26:41),2))./nanmean(rrr(any(rrr,2),26:41),2);
        imagesc(-0.2:0.01:0.3-0.01,[],rrrr)    
        colormap(brewermap([], '*RdBu'));
        caxis([-max(max(abs(rrrr))) max(max(abs(rrrr)))])
        
        figure(1);
        subplot(3, size(regions, 2), size(regions, 2)+(iRegion))
        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap);
        ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
        
        im = imagesc(xzCountBins{1, 1}./10, xzCountBins{1, 2}./10, binnedArrayPixelSmooth'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        %colorbar
        clearvars binnedArrayPixel
        bb2 = boundary(r, v, 0);
        hold on;
        plot(r(bb2)./10-(bregma(3) / 100), -v(bb2)./10, 'Color', theseColors{iRegion});
        xlabel('ML (mm)')
        ylabel('DV (mm)')
        axis equal
        axis square
        axis image
        
<<<<<<< HEAD
%         nColors = numel(ax.YTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
%         end
% 
%         nColors = numel(ax.XTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
%         end
%         ax.XLabel.Color = [0, 0, 0];
%         ax.YLabel.Color = [0, 0, 0];
=======
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
        ax.XLabel.Color = [0, 0, 0];
        ax.YLabel.Color = [0, 0, 0];
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        makepretty;
        clearvars isIN
        set(gca, 'color', [0.5, 0.5, 0.5])
        caxis(thisCmap)
        xlim([xzCountBins{1, 1}(1) ./ 10, xzCountBins{1, 1}(end) ./ 10])
        ylim([xzCountBins{1, 2}(1) ./ 10, xzCountBins{1, 2}(end) ./ 10])

        [N, Yedges, Zedges, binY, binZ] = histcounts2((theseLocations(:, 1)-bregma(1))./100, (-theseLocations(:, 2))./100, yzCountBins{1, 1}./10, yzCountBins{1, 2}./10); %par rapport a bregma!

        for iBinY = 1:size(Yedges, 2)
            for iBinZ = 1:size(Zedges, 2)
                theseNeurons = binY == iBinY & binZ == iBinZ;
                theseNeuronsInd = find(theseNeurons);
                binnedArrayTot = [];
                if ~isempty(theseNeuronsInd)
                    [h, hh] = histc(theseNeuronsInd, theseLocationsInfo);
                    uniqueRecs = unique(hh) + 1;
                    for iUniqueRec = 1:size(uniqueRecs, 1) %get psth per rec
                        theseTheseNeurons = theseNeuronsInd(hh == uniqueRecs(iUniqueRec));
                        if ~isempty(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(thisData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(thisData(uniqueRecs(iUniqueRec)).spike_templates, ...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec))-1);
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(thisData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                thisData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                            binnedArrayTot = [binnedArrayTot; binnedArray];
                            %figure();
                            %                    plot(-0.2:0.01:0.3-0.01,psth)
                        end
                    end

                end
                if ~isempty(binnedArrayTot)
                    binnedArrayTotBinned(iBinY, iBinZ, :) = nanmean(binnedArrayTot, 1); %average
                    % figure();
                    % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinY, iBinZ, :)))
                    if all(nanmean(binnedArrayTotBinned(iBinY, iBinZ, :)))
                        % figure();
                        % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinZ, :)))
                        binnedArrayPixel(iBinY, iBinZ) = (nanmean(binnedArrayTotBinned(iBinY, iBinZ, 26:41)) - ...
                            nanmean(binnedArrayTotBinned(iBinY, iBinZ, 1:20))) ./ nanmean(binnedArrayTotBinned(iBinY, iBinZ, 1:20));
                    else
                        binnedArrayTotBinned(iBinY, iBinZ, :) = nan(50, 1); %average
                        binnedArrayPixel(iBinY, iBinZ) = NaN;
                    end
                else
                    binnedArrayTotBinned(iBinY, iBinZ, :) = nan(50, 1); %average
                    binnedArrayPixel(iBinY, iBinZ) = NaN;
                end

            end
        end
<<<<<<< HEAD
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, smoothF, smoothF);
       % binnedArrayPixelSmooth = fliplr(binnedArrayPixelSmooth);
=======
        binnedArrayPixelSmooth = smooth2a(binnedArrayPixel, 4, 4);
        binnedArrayPixelSmooth = fliplr(binnedArrayPixelSmooth);
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
            for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
                isIN(iPixelX, iPixelY) = inpolygon(yzCountBins{1, 1}(iPixelX), ...
                    yzCountBins{1, 2}(iPixelY), c(bb3)-(bregma(1) / 10), -v(bb3));
            end
        end
        
        figure(iRegion+1)
        subplot(313)
        %plot(squeeze(binnedArrayTotBinned(iBinX, iBinY,:)))
        rrr = reshape(binnedArrayTotBinned(:, :,:), size(binnedArrayTotBinned,1).*size(binnedArrayTotBinned,2),size(binnedArrayTotBinned,3));
        rrrr = (rrr(any(rrr,2),:)-nanmean(rrr(any(rrr,2),26:41),2))./nanmean(rrr(any(rrr,2),26:41),2);
        imagesc(-0.2:0.01:0.3-0.01,[],rrrr)    
        colormap(brewermap([], '*RdBu'));
        caxis([-max(max(abs(rrrr))) max(max(abs(rrrr)))])
        
        figure(1)
        
        binnedArrayPixelSmooth(isIN == 0) = mean(thisCmap);
        subplot(3, size(regions, 2), size(regions, 2)*2+iRegion)
                ax = gca;
        ax.YColor = 'w'; % Red
        ax.XColor = 'w'; % Red
        
        im = imagesc(yzCountBins{1, 1}./10, yzCountBins{1, 2}./10, binnedArrayPixelSmooth'*100);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        clearvars isIN
<<<<<<< HEAD

        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
=======
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3

        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));

<<<<<<< HEAD
=======

>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
        %colorbar
        clearvars binnedArrayPixel
        subplot(3, size(regions, 2), size(regions, 2)*2+(iRegion))
        hold on;
        bb3 = boundary(c, v, 0);
        plot(c(bb3)./10-(bregma(1) / 100), -v(bb3)./10, 'Color', theseColors{iRegion});

        xlabel('AP (mm)')
        ylabel('DV (mm)')

        axis equal
        axis square
        axis image
        set(gca, 'color', [0.5, 0.5, 0.5]);

        %xticks([round(xl(1)):(round(xl(2))-round(xl(1))):round(xl(2))])
        %xticklabels({'x = 0','x = 5','x = 10'})
        makepretty;
        caxis(thisCmap)
        xlim([yzCountBins{1, 1}(1) ./ 10, yzCountBins{1, 1}(end) ./ 10])
        ylim([yzCountBins{1, 2}(1) ./ 10, yzCountBins{1, 2}(end) ./ 10])

<<<<<<< HEAD
%         nColors = numel(ax.YTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
%         end
%         ax.XLabel.Color = [0, 0, 0];
%         ax.YLabel.Color = [0, 0, 0];
%         nColors = numel(ax.XTickLabel);
%         cm = [0, 0, 0];
%         for i = 1:nColors
%             ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
%         end
=======
        nColors = numel(ax.YTickLabel);
        cm = [0, 0, 0];
        for i = 1:nColors
            ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
        end
        ax.XLabel.Color = [0, 0, 0];
        ax.YLabel.Color = [0, 0, 0];
        nColors = numel(ax.XTickLabel);
        cm = [0, 0, 0];
        for i = 1:nColors
            ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
        end
>>>>>>> a4b9a243a83c58bf74cc0891f6913b2d82bea8b3
set(gca, 'color', [0.5, 0.5, 0.5]);
    end
end

figure();
%binnedArrayPixelSmooth(isIN==0) = mean(thisCmap);
%   subplot(3, size(regions, 2), size(regions, 2)*2+iRegion)
im = imagesc(yzCountBins{1, 1}./10, yzCountBins{1, 2}./10, binnedArrayPixelSmooth'*100);
set(im, 'AlphaData', ~isnan(get(im, 'CData')));
clearvars isIN
colormap(brewermap([], '*RdBu'));
h = colorbar;
caxis(thisCmap)
h.Limits = [-5, 30];