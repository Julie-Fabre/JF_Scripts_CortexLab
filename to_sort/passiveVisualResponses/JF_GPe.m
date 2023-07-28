animal = 'AP083';
histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
load(histoFile)
probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
load(probe2ephysFile)
probe2ephys(5) = [];
recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
recInfo = readtable(recInfoFile);
allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);


protocolStrings = {'rating', 'ocation', 'nat'};
locations = {'GPi', 'PAL'};
%load probe
uniqueD = unique(recInfo.Date);
uniqueD = uniqueD(contains(uniqueD, '/')); %keep only dates (ie have a /) QQ
[dn, idx] = sort(datenum(uniqueD, 'dd/mm/yyyy'), 1, 'ascend'); %sort by ascending order
uniqueD = uniqueD(idx);

ephysData = struct;
thisCount = 1;

for iProbe = 3:4 %probe

    %get probe day and site.
    curr_day = probe2ephys(iProbe).day;
    curr_site = probe2ephys(iProbe).site;
    day = uniqueD{curr_day};
    dayBreaks = find(day == '/');
    if dayBreaks(1) == 2 && dayBreaks(2) == 4
        day = strcat(day(dayBreaks(2)+1:end), '-0', day(dayBreaks(1)+1:dayBreaks(2)-1), ...
            '-0', day(1:dayBreaks(1)-1)); %in standard format
    elseif dayBreaks(1) == 2
        day = strcat(day(dayBreaks(2)+1:end), '-', day(dayBreaks(1)+1:dayBreaks(2)-1), ...
            '-0', day(1:dayBreaks(1)-1)); %in standard format
    elseif dayBreaks(2) == 5
        day = strcat(day(dayBreaks(2)+1:end), '-0', day(dayBreaks(1)+1:dayBreaks(2)-1), ...
            '-', day(1:dayBreaks(1)-1)); %in standard format
    else
        day = strcat(day(dayBreaks(2)+1:end), '-', day(dayBreaks(1)+1:dayBreaks(2)-1), ...
            '-', day(1:dayBreaks(1)-1)); %in standard format
    end

    theseExpDates = strcmp(recInfo.Date, uniqueD{curr_day});
    if ~iscell(recInfo.Site)
        theseExpSites = recInfo.Site == curr_site;
    else
        a = recInfo.Site;
        b = cellfun(@numel, recInfo.Site);
        a(b > 1 | b == 0) = {'0'};
        c = str2num(cell2mat(a));
        theseExpSites = c == curr_site;
    end
    theseExpErrors = contains(recInfo.error, 'es') | ~isnan(recInfo.ephys_not_started) ...
        | ~isnan(recInfo.timeline_not_started);
    theseExp = find(theseExpDates & theseExpSites & ~theseExpErrors);
    theseProtocols = recInfo.Protocol(theseExp);


    for iProtocol = 1:3 %for each protocol type
        prot = theseExp(find(contains(theseProtocols, protocolStrings(iProtocol))));
        if ~isempty(prot)
            experiment = recInfo.Protocol_number(prot);
            site = c(prot);
            if numel(experiment) > 1
                prot = prot(end);
                experiment = experiment(end); %if several experiments, take the last one
                site = site(end);
            end


            AP_load_experimentJF;
            %get correct templates (location of interest)
            this_ccf = probe_ccf(iProbe);
            theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
            theseLocationsInterest = contains(theseLocations, locations{1}) | ...
                contains(theseLocations, locations{2});
            theseDepths = this_ccf.probe_depths(theseLocationsInterest);
            if ~isempty(theseDepths)
                %save data in structure
                theseTemplates = find(template_depths >= min(theseDepths) & template_depths <= max(theseDepths)); %correct depth units
                %find closest depth and get closest location
                A = repmat(this_ccf.probe_depths, [1, length(template_depths)]);
                [minValue, closestIndex] = min(abs(A-template_depths'));
                closestLocation = this_ccf.trajectory_coords(closestIndex, :);
                %raster general and for each trial condition
                theseSpikes = ismember(spike_templates, theseTemplates);
                %save info
                ephysData(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
                ephysData(thisCount).spike_templates = spike_templates(theseSpikes);
                ephysData(thisCount).templatesOI = theseTemplates;
                ephysData(thisCount).template_depths = template_depths(theseTemplates);
                ephysData(thisCount).template_wf = waveforms(theseTemplates, :);
                ephysData(thisCount).template_location = closestLocation;
                ephysData(thisCount).protocol = expDef;
                ephysData(thisCount).animal = animal;
                ephysData(thisCount).date = day;
                ephysData(thisCount).site = iProbe;
                ephysData(thisCount).experiment = experiment;
                ephysData(thisCount).location = 'GPi';
                ephysData(thisCount).stimOn_times = stimOn_times;
                ephysData(thisCount).stimIDs = stimIDs;
                ephysData(thisCount).trial_conditions = trial_conditions;
                thisCount = thisCount + 1;
            end
        end
    end
end
%% averages -all

    colorsT = {'r', 'r', 'r', 'g', 'g', 'g', 'b', 'b', 'b'};
    figure();
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    for i = [1, 4]
        theseSpikes = ephysData(i).spike_times_timeline;
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, [-0.2, 0], psthBinSize);
        psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT{i}, 0.1)
        makepretty;
        hold on;
    end

    ylabel('dFR/FR')
    xlabel('time from grating (s)')
    set(gcf, 'color', 'w');

        colorsT = {' ',  ' ', 'r', 'r', ' ', 'g'};
    figure();
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    for i = [3, 6]
        theseSpikes = ephysData(i).spike_times_timeline;
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, [-0.2, 0], psthBinSize);
        psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT{i}, 0.1)
        makepretty;
        hold on;
    end

    ylabel('dFR/FR')
    xlabel('time from nat. image (s)')
    set(gcf, 'color', 'w');
%% averages - orientations
% trial_conditions = ...
%                [signals_events.stimAzimuthValues', signals_events.stimSpatialFreqValues', ...
%                signals_events.stimOrientationValues'];

%psth for each rec to location, nat image, grating
for i = [1, 4]
uniq = unique(ephysData(1).trial_conditions(:,3));
    colorsT = {rgb('Purple'), rgb('DarkRed'), rgb('Red'), rgb('OrangeRed'), rgb('Orange'), rgb('Gold'), rgb('YellowGreen'), rgb('Green')};
    figure();
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    for iOrien = 1:8
        theseSpikes = ephysData(i).spike_times_timeline;
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).trial_conditions(:,3)==uniq(iOrien)), thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).trial_conditions(:,3)==uniq(iOrien)), [-0.2, 0], psthBinSize);
        psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT{iOrien},0.1)
        makepretty;
        hold on;
    end

    ylabel('dFR/FR')
    xlabel('time from grating (s)')
    set(gcf, 'color', 'w');
end

for i = [1, 4]
uniq = unique(ephysData(1).trial_conditions(:,2));
    colorsT = {rgb('Purple'), rgb('DarkRed'), rgb('Red'), rgb('OrangeRed'), rgb('Orange'), rgb('Gold'), rgb('YellowGreen'), rgb('Green')};
    figure();
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    for iOrien = 1:4
        theseSpikes = ephysData(i).spike_times_timeline;
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).trial_conditions(:,2)==uniq(iOrien)), thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).trial_conditions(:,2)==uniq(iOrien)), [-0.2, 0], psthBinSize);
        psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT{iOrien},0.1)
        makepretty;
        hold on;
    end

    ylabel('dFR/FR')
    xlabel('time from grating (s)')
    set(gcf, 'color', 'w');
end

for i = [2, 5]
    uniq = unique(ephysData(2).stimIDs);
    colorsT = {rgb('MediumVioletRed'); rgb('DeepPink'); rgb('HotPink'); rgb('Brown'); rgb('Purple');rgb('DarkBlue');rgb('Blue');rgb('SkyBlue')};
    figure();
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    for iOrien = 1:6
        theseSpikes = ephysData(i).spike_times_timeline;
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).stimIDs==uniq(iOrien)), thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(ephysData(i).stimIDs==uniq(iOrien)), [-0.2, 0], psthBinSize);
        psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT{iOrien},0.1)
        makepretty;
        hold on;
    end

    ylabel('dFR/FR')
    xlabel('time from grating (s)')
    set(gcf, 'color', 'w');
end

thisWindow = [-0.2, 0.5];
psthBinSize = 0.01;
ccount = 0;
for i = [3,6]

        oo =1:30;
        ooo = 1:length(oo);
    
    for iOr = ooo
        
            ff = find(ephysData(i).stimIDs == oo(iOr));
        
        
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(ephysData(i).spike_times_timeline, ephysData(i).stimOn_times(ff), thisWindow, psthBinSize);
        [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(ephysData(i).spike_times_timeline, ephysData(i).stimOn_times(ff), [-0.2, 0], psthBinSize);
        locationMUAPSTH(i, iOr, :) = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
        w = gausswin(10);
        y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
        locationMUAPSTH_STD(i, iOr, :) = std(y) ./ sqrt(size(binnedArray, 1));
    end
    figure();
imagesc([-0.2:0.01:0.5 - 0.01]', [],squeeze((locationMUAPSTH(i,  :,:))))
colormap(gray)

xlabel('time from stim onset (s)')
ylabel('image #')
yy=ylim;
line([0 0],[yy(1) yy(2)],'Color','r')
h = colorbar;
ylabel(h, 'dFR/FR')
%ylim([-4, 6])
makeprettyLite;
end

figure();
imagesc([-0.2:0.01:0.5 - 0.01]', [],squeeze(nanmean(locationMUAPSTH(1:5,  :,:))))
colormap(gray)

xlabel('time from stim onset (s)')
ylabel('image #')
%ylim([-4, 6])
makeprettyLite;

%% individual cells : locations

figure();
iRec = 5;
iTemplate = 1;

iTemplate=iTemplate+1
colorsO = [rgb('Red'); rgb('BlueViolet'); rgb('Purple'); rgb('Blue'); rgb('DeepSkyBlue'); rgb('Turquoise'); rgb('MediumSeaGreen'); rgb('Green')];
uu = unique(ephysData(iRec).spike_templates);
thisWindow = [-0.2, 0.5];
psthBinSize = 0.01;
iTemplate = iTemplate + 1;

thisT = uu(iTemplate);
clf;
binnFull = [];
iu = unique(ephysData(iRec).stimIDs);
for iiu = 1:length(iu)
    theseSpikes = ephysData(iRec).spike_times_timeline(ephysData(iRec).spike_templates == thisT);
    [psth, ~, ~, ~, ~, binnedArray] = psthAndBA(theseSpikes, ephysData(iRec).stimOn_times(ephysData(iRec).stimIDs == iu(iiu)), thisWindow, psthBinSize);
    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(iRec).stimOn_times(ephysData(iRec).stimIDs == iu(iiu)), [-0.2, 0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    subplot(4, 1, 4)
    w = gausswin(10);
    y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
    ss = std(y) ./ sqrt(size(binnedArray, 1));
    AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', ss, colorsO(iiu, :))
    makepretty;
    xlim([-0.2, 0.5])
    hold on;

    orTun(iiu) = nanmean(psth(:, 20:40));
    orTunSTD(iiu) = nanstd(psth(:, 20:40));
    binnFull = [binnFull; binnedArray];
end

subplot(4, 1, [1:3])
binnFull(binnFull > 1) = 1;
h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnFull);
colormap(gray)
yy = ylim;
for iOr = 1:length(iu)
    line([0, 0], [1 + round(yy(2)/length(iu)) * (iOr - 1), round(yy(2)/length(iu)) * (iOr)], 'Color', colorsO(iOr, :))
    makepretty;
end
disp(iTemplate)

%% inidividual cells: orientation

%% inidividual cells: location 

%% individual cells: anything??
figure();
iRec = 2;
iTemplate = 1;

iTemplate=iTemplate+1
colorsO = [rgb('Red'); rgb('BlueViolet'); rgb('Purple'); rgb('Blue'); rgb('DeepSkyBlue'); rgb('Turquoise'); rgb('MediumSeaGreen'); rgb('Green')];
uu = unique(ephysData(iRec).spike_templates);
thisWindow = [-0.2, 0.5];
psthBinSize = 0.01;


thisT = uu(iTemplate);
clf;
binnFull = [];
iu = unique(ephysData(iRec).stimIDs);
 theseSpikes = ephysData(iRec).spike_times_timeline(ephysData(iRec).spike_templates == thisT);
    [psth, ~, ~, ~, ~, binnedArray] = psthAndBA(theseSpikes, ephysData(iRec).stimOn_times, thisWindow, psthBinSize);
    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(iRec).stimOn_times, [-0.2, 0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    subplot(4, 1, 4)
    w = gausswin(10);
    y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
    ss = std(y) ./ sqrt(size(binnedArray, 1));
    AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', ss, colorsO(iiu, :))
    makepretty;
    xlim([-0.2, 0.5])
    hold on;
     orTun(iiu) = nanmean(psth(:, 20:40));
    orTunSTD(iiu) = nanstd(psth(:, 20:40));
    binnFull = [binnFull; binnedArray];
    makeprettyLite; 
    
subplot(4, 1, [1:3])
binnFull(binnFull > 1) = 1;
h = imagesc(-0.2:0.01:0.5-0.01, [], 1-binnFull);
colormap(gray)
yy = ylim;

%for iiu = 1:length(iu)
   

   
%end

disp(iTemplate)