
%% load naive visual recordings
corona = 0;
mice = {'AP080', 'AP082', 'JF017', 'JF019'}; % JF020 doesn't have histology yet. JF021 and JF022 to come, and many more
locations = {'CP', 'SNr', 'GPe', 'GPi', 'STN'};
protocols = {'JF_natural_images', 'JF_locations', 'JF_GratingsPassive'};
ephysData = struct;

allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);
thisCount = 1;
for iMouse = 1:size(mice, 2)
    %find experiments
    animal = mice{iMouse};
    protocol = 'JF'; %contains this - refine protocol later
    experiments = AP_find_experimentsJF(animal, protocol, 1);
    experiments = experiments([experiments.ephys]);
    % check if kilosorted
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)

    recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
    recInfo = readtable(recInfoFile);
    for iDate = 1:size(experiments, 1)
        corrDay = find(arrayfun(@(x) probe2ephys(x).day == iDate, 1:numel(probe2ephys)));
        for iCorrDay = 1:size(corrDay, 2)
            sites(iCorrDay) = probe2ephys(corrDay(iCorrDay)).site;
        end
        for iSite = 1:size(sites, 2)
            corrSite = find(arrayfun(@(x) probe2ephys(x).site == iSite, 1:numel(probe2ephys)));
            for iLocation = 1:size(locations, 2)
                probe = intersect(corrDay, corrSite);
                %disp(probe)
                thisP = probe;
                if ~isempty(probe)
                    this_ccf = probe_ccf(probe); %find if contains location of interest
                    if ~any(structfun(@isempty, this_ccf)) %check hewre
                        theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
                        theseLocationsInterest = contains(theseLocations, locations{iLocation});
                        theseDepths = this_ccf.probe_depths(theseLocationsInterest);
                        if ~isempty(find(theseDepths)) %check hewre
                            %load experiment
                            curr_day = probe2ephys(probe).day;
                            day = experiments(curr_day).day;
                            experimentThese = experiments(curr_day).experiment; % experiment number

                            % remove experiments where bug (eg
                            % timeline/ephys not started)
                            thisProtocol = find(ismember(recInfo.Protocol_number, [experimentThese]) & ...
                                ismember(recInfo.Date, day));
                            %
                            if iscell(recInfo.ephys_not_started(1))
                                ephysOK = cellfun(@isempty, recInfo.ephys_not_started(thisProtocol));
                            else
                                ephysOK = ones(length(experimentThese), 1);
                            end
                            if iscell(recInfo.timeline_not_started(1))
                                timelineOK = cellfun(@isempty, recInfo.timeline_not_started(thisProtocol));
                            else
                                timelineOK = ones(length(experimentThese), 1);
                            end

                            experimentThese = experimentThese(ephysOK & timelineOK);


                            %experiment = experiment(probe2ephys(probe).site);
                            verbose = false; % display load progress and some info figures
                            load_parts.cam = false;
                            load_parts.imaging = false;
                            load_parts.ephys = true;
                            site = probe2ephys(probe).site;
                            for iExperiment = 1:size(experimentThese, 2)
                                experiment = experimentThese(iExperiment);
                                AP_load_experimentJF;
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
                                try
                                    ephysData(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
                                end
                                ephysData(thisCount).templatesOI = theseTemplates;
                                ephysData(thisCount).template_depths = template_depths(theseTemplates);

                                ephysData(thisCount).template_wf = waveforms(theseTemplates, :);

                                ephysData(thisCount).template_location = closestLocation;
                                ephysData(thisCount).protocol = expDef;
                                ephysData(thisCount).animal = animal;
                                ephysData(thisCount).date = day;
                                ephysData(thisCount).site = thisP;
                                ephysData(thisCount).experiment = experiment;
                                ephysData(thisCount).location = locations(iLocation);
                                ephysData(thisCount).stimOn_times = stimOn_times;
                                ephysData(thisCount).stimIDs = stimIDs;
                                thisCount = thisCount + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% plot responses - Striatum

%find == 'CP'
striatumData = find(arrayfun(@(x) contains(ephysData(x).location, 'CP'), 1:numel(ephysData)));

%find striatum limits - allen atlas
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure}, 2, [])') ./ 255;

%get limits and bins
binSize = 100;
structure_3d_lims = permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]);
[r, c, v] = ind2sub(size(structure_3d_lims), find(structure_3d_lims));
x_lim_um = [min(r), max(r)] .* 100;
y_lim_um = [min(c), max(c)] .* 100;
z_lim_um = [min(v), max(v)] .* 100;
dataBinsX = x_lim_um(1):binSize:x_lim_um(2);
dataBinsY = y_lim_um(1):binSize:y_lim_um(2);
dataBinsZ = z_lim_um(1):binSize:z_lim_um(2);

%find which bins striatal templates belong to
uniqueRecordings = arrayfun(@(x) [ephysData(x).animal, ephysData(x).date, num2str(ephysData(x).site)], ...
    1:numel(ephysData), 'UniformOutput', false); %unique date + site
[uniqueId, uniqueIdx] = unique(uniqueRecordings, 'rows');
%psth - combining all protocols
posBinSize = 50;
timeBinSize = 0.01;
win = [0, 0.5];
bslWin = [-0.2, 0];

allP_allprotocols_allrecs = [];
posBinsX_allprotocols_allrecs = [];
posBinsY_allprotocols_allrecs = [];

for iUniqueRec = 1:size(uniqueId, 2)
    theseRecs = find(contains(uniqueRecordings, uniqueId(iUniqueRec)));
    allP_allprotocols = [];
    posBinsX_allprotocols = [];
    posBinsY_allprotocols = [];
    for iProtocol = 1:size(theseRecs, 2)

        spikeTimes = ephysData(theseRecs(iProtocol)).spike_times_timeline;
        uniqueTemp = unique(ephysData(theseRecs(iProtocol)).spike_templates);
        spikePos = nan(size(spikeTimes, 1), 3);
        for iTemp = 1:size(uniqueTemp, 1)
            theseSp = find(ephysData(theseRecs(iProtocol)).spike_templates == uniqueTemp(iTemp));
            spikePos(theseSp, 1) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 1);
            spikePos(theseSp, 2) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 2);
            spikePos(theseSp, 3) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 3);
        end
        eventTimes = ephysData(theseRecs(iProtocol)).stimOn_times;

        if (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize && (max(spikePos(:, 2)) - min(spikePos(:, 2))) > posBinSize
            [timeBins, posBinsX, posBinsY, allP, normVals] = psthByPos2D(spikeTimes, spikePos(:, 1), squeeze(spikePos(:, 2)), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);

        elseif (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize
            [timeBins, posBinsX, allP, normVals] = psthByPos1D(spikeTimes, spikePos(:, 1), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsY = nan(size(posBinsX, 2), 1);
            posBinsY(:) = nanmean(spikePos(:, 2));

        elseif (max(spikePos(:, 2)) - min(spikePos(:, 2))) > posBinSize
            [timeBins, posBinsY, allP, normVals] = psthByPos1D(spikeTimes, squeeze(spikePos(:, 2)), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsX = nan(size(posBinsY, 2), 1);
            posBinsX(:) = nanmean(spikePos(:, 1));
        else
            [psth, ~, ~, ~, ~, ~] = psthAndBA(spikeTimes, eventTimes, bslWin, timeBinSize);
            normMn = mean(psth);
            normStd = std(psth);

            [psth, ~, ~, ~, ~, ~] = psthandBA(spikeTimes, ...
                eventTimes, win, timeBinSize);
            allP = (psth - normMn) ./ normStd;
            posBinsX = nanmean(spikePos(:, 1));
            posBinsY = nanmean(spikePos(:, 2));
        end
        allP_allprotocols = [allP_allprotocols, allP];
        posBinsX_allprotocols = [posBinsX_allprotocols, posBinsX];
        posBinsY_allprotocols = [posBinsY_allprotocols, posBinsY];
    end
    allP_allprotocols_allrecs = [allP_allprotocols_allrecs, allP_allprotocols];
    posBinsX_allprotocols_allrecs = [posBinsX_allprotocols_allrecs, posBinsX_allprotocols];
    posBinsY_allprotocols_allrecs = [posBinsY_allprotocols_allrecs, posBinsY_allprotocols];

end

%psth per grating, location, orientation, sp. frequency, nat. image

%cell type
