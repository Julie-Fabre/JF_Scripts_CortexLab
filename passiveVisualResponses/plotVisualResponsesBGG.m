keep ephysData

animalsType = {'Naive'};
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
%regions = {'DMS', 'PS', 'STN', 'GPe', 'SNr', 'GPi'};
recordingInfo = readtable('C:\Users\Julie\Dropbox\Analysis\Recordings - Sheet1.csv');
bregma = [540, 0, 570];

%% X-Y, X-Z, Y-Z plots for each region. bin by 100 um (?) and plot in colormap increase in response ?

for iType = 1:size(animalsType, 2)
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
    for iRegion = 1:size(regions, 2)

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
        bb = boundary(r, c);
        plot(r(bb), c(bb), 'Color', theseColors{iRegion});
        xlabel('ML')
        ylabel('AP')
        axis equal
        axis square
        xyXLim = xlim;
        xyYLim = ylim;
        makepretty;

        % x-z
        subplot(3, size(regions, 2), size(regions, 2)+iRegion)
        bb = boundary(r, v);
        plot(r(bb), v(bb), 'Color', theseColors{iRegion});
        xlabel('ML')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
         xzXLim = xlim;
        xzZLim = ylim;

        % y-z
        subplot(3, size(regions, 2), size(regions, 2)*2+(iRegion))
        bb = boundary(c, v);
        plot(c(bb) , v(bb), 'Color', theseColors{iRegion});
        xlabel('AP')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
        yzYLim = xlim;
        yzZLim = ylim;
        % bin region in 200um bins
        xy = [r(bb), c(bb)];
        xyCountBins = {xyXLim(1) * 10 :10:xyXLim(2) * 10 , xyYLim(1) * 10 :10:xyYLim(2) * 10 };
        xzCountBins = {xzXLim(1) * 10 :10:xzXLim(2) * 10 , xzZLim(1) * 10 :10:xzZLim(2) * 10 };
        yzCountBins = {yzYLim(1) * 10 :10:yzYLim(2) * 10 , yzZLim(1) * 10 :10:yzZLim(2) * 10 };
        % load all region data [-0.1 - 0.3]s after stim onset, with which bin belongs to

        CNew = cat(1, ephysData.location);

        ind = find(ismember(CNew, regions{iRegion}));

        for iLocation = 1:size(ind, 1)
            theseLocationsInfo(iLocation) = size(ephysData(ind(iLocation)).template_location, 1);
        end
        theseLocationsInfo = cumsum(theseLocationsInfo);
        theseLocations = cat(1, ephysData(ind).template_location);
        theseLocationsBregmaAbs = [abs(theseLocations(:, 3)-bregma(3)), abs(theseLocations(:, 1)-bregma(1)), theseLocations(:, 2)];
        
        subplot(3, size(regions, 2), iRegion)
        hold on;
        scatter((bregma(3)/10) -abs(theseLocations(:,3)/10 - (bregma(3)/10)), theseLocations(:,1)/10)
        
%         figure();
%         hold on;
%         scatter(theseLocations(:,3)/10, theseLocations(:,1)/10)
%         scatter((bregma(3)/10) -abs(theseLocations(:,3)/10 - (bregma(3)/10)), theseLocations(:,1)/10)
%         scatter(bregma(1)/10, bregma(3)/10)
%         scatter(bregma(3)/10, bregma(1)/10)
%         hold on;
%         
        % x-z
        subplot(3, size(regions, 2), size(regions, 2)+iRegion)
        hold on;
        scatter((bregma(3)/10) -abs(theseLocations(:,3)/10 - (bregma(3)/10)), theseLocations(:,2)/10)

        % y-z
        subplot(3, size(regions, 2), size(regions, 2)*2+(iRegion))
        hold on;
        scatter(theseLocations(:,1)/10, theseLocations(:,2)/10)
         subplot(3, size(regions, 2), iRegion)
        bb = boundary(r, c);
        plot(r(bb), c(bb), 'Color', theseColors{iRegion});
        xlabel('ML')
        ylabel('AP')
        axis equal
        axis square
        xyXLim = xlim;
        xyYLim = ylim;
        makepretty;

        % x-z
        subplot(3, size(regions, 2), size(regions, 2)+iRegion)
        bb = boundary(r, v);
        plot(r(bb), v(bb), 'Color', theseColors{iRegion});
        xlabel('ML')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
         xzXLim = xlim;
        xzZLim = ylim;

        % y-z
        subplot(3, size(regions, 2), size(regions, 2)*2+(iRegion))
        bb = boundary(c, v);
        plot(c(bb) , v(bb), 'Color', theseColors{iRegion});
        xlabel('AP')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
        yzYLim = xlim;
        yzZLim = ylim;
        
        %
%        theseSpikeTimes = cat(1, ephysData(ind).spike_times_timeline);
%        theseSpikeTemplates = cat(1, ephysData(ind).spike_templates);
%        theseStimOnTimes = cat(ephysData(ind).stimOn_times, 1);
%         hold on;
%        plot3(theseLocations(:, 1), theseLocations(:, 3), theseLocations(:, 2), '.b') %1,2,3;23,2,1;
        [N, Xedges, Yedges, binX, binY] = histcounts2((bregma(3)) -abs(theseLocations(:,3) - (bregma(3))), theseLocations(:, 1), xyCountBins{1, 1}, xyCountBins{1, 2}); %par rapport a bregma!

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
                        if ~isempty(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(ephysData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(ephysData(uniqueRecs(iUniqueRec)).spike_templates, ...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec)));
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                ephysData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                            binnedArrayTot = [binnedArrayTot; binnedArray];
%figure();
%                    plot(-0.2:0.01:0.3-0.01,psth)
                        end
                    end

                end
                if ~isempty(binnedArrayTot)
                    binnedArrayTotBinned(iBinX, iBinY, :) = nanmean(binnedArrayTot, 1); %average
                    % figure(); 
                    % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinY, :)))
                    binnedArrayPixel(iBinX, iBinY) = (nanmean(binnedArrayTotBinned(iBinX, iBinY, 26:41)) - ...
                        nanmean(binnedArrayTotBinned(iBinX, iBinY, 1:20))) ./nanmean(binnedArrayTotBinned(iBinX, iBinY, 1:20));
                else
                    binnedArrayTotBinned(iBinX, iBinY, :) = nan(50, 1); %average
                end
            end
        end
        figure(1+iRegion);
        subplot(131)
        imagesc(binnedArrayPixel*100)
        colorbar
        clearvars binnedArrayPixel 
        
        [N, Xedges, Zedges, binX, binZ] = histcounts2((bregma(3)) -abs(theseLocations(:,3) - (bregma(3))), theseLocations(:, 2), xzCountBins{1, 1}, xzCountBins{1, 2}); %par rapport a bregma!

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
                        if ~isempty(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(ephysData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(ephysData(uniqueRecs(iUniqueRec)).spike_templates, ...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec)));
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                ephysData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                            binnedArrayTot = [binnedArrayTot; binnedArray];
%figure();
%                    plot(-0.2:0.01:0.3-0.01,psth)
                        end
                    end

                end
                if ~isempty(binnedArrayTot)
                    binnedArrayTotBinned(iBinX, iBinZ, :) = nanmean(binnedArrayTot, 1); %average
                    % figure(); 
                    % plot(-0.2:0.01:0.3-0.01,squeeze(binnedArrayTotBinned(iBinX, iBinZ, :)))
                    binnedArrayPixel(iBinX, iBinZ) = (nanmean(binnedArrayTotBinned(iBinX, iBinZ, 26:41)) - ...
                        nanmean(binnedArrayTotBinned(iBinX, iBinZ, 1:20))) ./nanmean(binnedArrayTotBinned(iBinX, iBinZ, 1:20));
                else
                    binnedArrayTotBinned(iBinX, iBinZ, :) = nan(50, 1); %average
                end
            end
          end
          figure(1+iRegion);
        subplot(132)
        imagesc(binnedArrayPixel*100)
        colorbar
        clearvars binnedArrayPixel 
            [N, Yedges, Zedges, binY, binZ] = histcounts2(theseLocations(:, 1), theseLocations(:, 2), yzCountBins{1, 1}, yzCountBins{1, 2}); %par rapport a bregma!

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
                        if ~isempty(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(ephysData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(ephysData(uniqueRecs(iUniqueRec)).spike_templates, ...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec)));
                            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                psthAndBA(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate), ...
                                ephysData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
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
                    binnedArrayPixel(iBinY, iBinZ) = (nanmean(binnedArrayTotBinned(iBinY, iBinZ, 26:41)) - ...
                        nanmean(binnedArrayTotBinned(iBinY, iBinZ, 1:20))) ./nanmean(binnedArrayTotBinned(iBinY, iBinZ, 1:20));
                else
                    binnedArrayTotBinned(iBinY, iBinZ, :) = nan(50, 1); %average
                end
            end
          end
          figure(1+iRegion);
        subplot(133)
        imagesc(binnedArrayPixel*100)
        colorbar
        clearvars binnedArrayPixel 
    end
end