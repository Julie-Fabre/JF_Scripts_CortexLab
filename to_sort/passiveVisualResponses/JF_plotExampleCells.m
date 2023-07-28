figure(8);
clf;

iRegion = 1;
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
thisData = ephysData;

iRegion = 4;
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
thisData = ephysData;

%iUnit = 1;
%all Units
CNew = cat(1, thisData.location);
CProt = strcat(thisData.protocol);
newStr = split(CProt, 'J');
load([dropboxPath 'Presentations/UpgradeDataAndPlots/conditionsOriGratPassive.mat'])
load([dropboxPath 'Presentations/UpgradeDataAndPlots/conditionsOriGratPassiveNew.mat'])

ind = find(ismember(CNew, regions{iRegion}) & contains(newStr(2:end), 'ratingPassive'));


thisRec = thisRec - 1;

thisRec=13
iUnit = 0;
theseTemplates = unique(thisData(ind(thisRec)).spike_templates);

%% Loop
clf;
iUnit = iUnit + 1;
theseTemplates = unique(thisData(ind(thisRec)).spike_templates);

if contains(newStr(ind(thisRec)+2), 'natural_image') %previous, locations in gratings1 scripr
    allOri = unique(conditions(:, 3));
    binnedArrayFull = [];
    binnedArrayType = [];
    for iOri = 1:size(allOri, 1)
        allOriRows = conditions(:, 3) == allOri(iOri);
        [psthOr(iOri, :), bins, rasterX, rasterY, spikeCounts, binnedArrayFull] = ...
            psthAndBA(thisData(ind(thisRec)).spike_times_timeline(thisData(ind(thisRec)).spike_templates == theseTemplates(iUnit)), ...
            thisData(ind(thisRec)).stimOn_times(thisData(ind(thisRec)).stimIDs == allOriRows), [-0.2, 0.5], 0.01);
        binnedArrayFull = [binnedArrayFull, binnedArrayFull];
        binnedArrayType = [binnedArrayType, ones(size(binnedArrayFull, 1), 1) .* iOri];
    end
    subplot(6, 2, [1:3])

    [psth, bins, rasterX, rasterY, spikeCounts, binnedArrayFull] = ...
        psthAndBA(thisData(ind(thisRec)).spike_times_timeline(thisData(ind(thisRec)).spike_templates == theseTemplates(iUnit)), ...
        thisData(ind(thisRec)).stimOn_times, [-0.2, 0.5], 0.01);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArrayFull] = ...
        psthAndBA(thisData(ind(thisRec)).spike_times_timeline(thisData(ind(thisRec)).spike_templates == theseTemplates(iUnit)), ...
        thisData(ind(thisRec)).stimOn_times, [-0.2, 0.5], 0.01);
else
    %% orientations
    theseTemplates = unique(thisData(ind(thisRec)).spike_templates);

    allOri = unique(conditionsNew(:, 3));
    binnedArrayFull = [];
    binnedArrayType = [];
    for iOri = 1:size(allOri, 1)
        allOriRows = find(conditionsNew(:, 3) == allOri(iOri));
        [psthOr(iOri, :), bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
            psthAndBA(thisData(ind(thisRec)).spike_times_timeline(thisData(ind(thisRec)).spike_templates == theseTemplates(iUnit)), ...
            thisData(ind(thisRec)).stimOn_times(ismember(thisData(ind(thisRec)).stimIDs, allOriRows)), [-0.2, 0.5], 0.01);
        binnedArrayFull = [binnedArrayFull, binnedArray'];
        binnedArrayType = [binnedArrayType, ones(size(binnedArrayFull, 1), 1) .* iOri];
    end
    subplot(6, 2, [1, 3])
     binnedArrayFull=binnedArrayFull';
    binnedArrayFull(binnedArrayFull > 1) = 1;
    h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnedArrayFull);
    colormap(gray)
    yy = ylim;
    colorsO = {rgb('DarkRed'), rgb('Red'), rgb('OrangeRed'), rgb('DarkOrange'), rgb('GoldenRod'), rgb('LimeGreen'), rgb('ForestGreen'), rgb('Teal'), rgb('Purple')};
    for iOri = 1:size(allOri, 1)
        line([0, 0], [round(yy(2)/size(allOri, 1)) .* (iOri - 1) + 1, yy(2)], 'Color', colorsO{iOri})
    end
    makepretty;

    subplot(6, 2, 5)
    for iOri = 1:size(allOri, 1)

        plot(-0.2:0.01:0.5-0.01, squeeze(psthOr(iOri, :)), 'Color', colorsO{iOri})
        hold on;
        makepretty;
        xlim([-0.2, 0.5])
    end
    %% spatial frequency 
     allOri = unique(conditionsNew(:, 2));
    binnedArrayFull = [];
    binnedArrayType = [];
    for iOri = 1:size(allOri, 1)
        allOriRows = find(conditionsNew(:, 2) == allOri(iOri));
        [psthOr(iOri, :), bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
            psthAndBA(thisData(ind(thisRec)).spike_times_timeline(thisData(ind(thisRec)).spike_templates == theseTemplates(iUnit)), ...
            thisData(ind(thisRec)).stimOn_times(ismember(thisData(ind(thisRec)).stimIDs, allOriRows)), [-0.2, 0.5], 0.01);
        binnedArrayFull = [binnedArrayFull, binnedArray'];
        binnedArrayType = [binnedArrayType, ones(size(binnedArrayFull, 1), 1) .* iOri];
    end
    subplot(6, 2, [2, 4])
     binnedArrayFull=binnedArrayFull';
    binnedArrayFull(binnedArrayFull > 1) = 1;
    h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnedArrayFull);
    colormap(gray)
    yy = ylim;
    colorsO = {rgb('DarkRed'), rgb('Red'), rgb('OrangeRed'), rgb('DarkOrange'), rgb('GoldenRod'), rgb('LimeGreen'), rgb('ForestGreen'), rgb('Teal'), rgb('Purple')};
    for iOri = 1:size(allOri, 1)
        line([0, 0], [round(yy(2)/size(allOri, 1)) .* (iOri - 1) + 1, yy(2)], 'Color', colorsO{iOri})
    end
    makepretty;

    subplot(6, 2, 6)
    for iOri = 1:size(allOri, 1)

        plot(-0.2:0.01:0.5-0.01, squeeze(psthOr(iOri, :)), 'Color', colorsO{iOri})
        hold on;
        makepretty;
        xlim([-0.2, 0.5])
    end
    %% location
    if contains(thisData(ind(thisRec)+1).protocol, 'ocation')
        addThis = 0;
    else
        addThis = 1;
    end
    theseTemplates = unique(thisData(ind(thisRec)+1+addThis).spike_templates);

    allOri = 1:6;
    binnedArrayFull = [];
    binnedArrayType = [];
    for iOri = 1:size(allOri, 2)
        allOriRows = iOri;
        [psthOr(iOri, :), bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
            psthAndBA(thisData(ind(thisRec)+1+addThis).spike_times_timeline(thisData(ind(thisRec)+1+addThis).spike_templates == theseTemplates(iUnit)), ...
            thisData(ind(thisRec)+1+addThis).stimOn_times(ismember(thisData(ind(thisRec)+1+addThis).stimIDs, allOriRows)), [-0.2, 0.5], 0.01);
        binnedArrayFull = [binnedArrayFull, binnedArray'];
        binnedArrayType = [binnedArrayType, ones(size(binnedArrayFull, 1), 1) .* iOri];
    end
    subplot(6, 2, [7, 9])
     binnedArrayFull=binnedArrayFull';
    binnedArrayFull(binnedArrayFull > 1) = 1;
    h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnedArrayFull);
    colormap(gray)
    yy = ylim;
    colorsO = {rgb('DarkRed'), rgb('Red'), rgb('OrangeRed'), rgb('DarkOrange'), rgb('GoldenRod'), rgb('LimeGreen'), rgb('ForestGreen'), rgb('Teal'), rgb('Purple')};
    for iOri = 1:size(allOri, 2)
        line([0, 0], [round(yy(2)/size(allOri, 2)) .* (iOri - 1) + 1, yy(2)], 'Color', colorsO{iOri})
    end
    makepretty;

    subplot(6, 2, 11)
    for iOri = 1:size(allOri, 2)

        plot(-0.2:0.01:0.5-0.01, squeeze(psthOr(iOri, :)), 'Color', colorsO{iOri})
        hold on;
        makepretty;
        xlim([-0.2, 0.5])
    end
    %% natural images 
     if contains(thisData(ind(thisRec)+2+addThis).protocol, 'atural_images')
        addThis = addThis;
    else
        addThis = addThis+1;
     end
    
    theseTemplates = unique(thisData(ind(thisRec)+2+addThis).spike_templates);

    allOri = 1:30;
    binnedArrayFull = [];
    binnedArrayType = [];
    for iOri = 1:size(allOri, 2)
        allOriRows = iOri;
        [psthOr(iOri, :), bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
            psthAndBA(thisData(ind(thisRec)+2+addThis).spike_times_timeline(thisData(ind(thisRec)+2+addThis).spike_templates == theseTemplates(iUnit)), ...
            thisData(ind(thisRec)+2+addThis).stimOn_times(ismember(thisData(ind(thisRec)+2+addThis).stimIDs, allOriRows)), [-0.2, 0.5], 0.01);
        binnedArrayFull = [binnedArrayFull, binnedArray'];
        binnedArrayType = [binnedArrayType, ones(size(binnedArrayFull, 1), 1) .* iOri];
    end
    subplot(6, 2, [8, 10])
    binnedArrayFull=binnedArrayFull';
    binnedArrayFull(binnedArrayFull > 1) = 1;
    h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnedArrayFull);
    colormap(gray)
    yy = ylim;
    colorsO = lines(30);
    for iOri = 1:size(allOri, 2)
        line([0, 0], [round(yy(2)/size(allOri, 2)) .* (iOri - 1) + 1, yy(2)], 'Color', colorsO(iOri,:))
    end
    makepretty;

    subplot(6, 2, 12)
    imagesc(-0.2:0.01:0.5-0.01, [], psthOr)
    makepretty;
end

disp(thisRec)
disp(iUnit)
