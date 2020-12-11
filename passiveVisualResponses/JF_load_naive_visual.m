
%% load naive visual recordings
corona = 0;
mice = {'AP080', 'AP082', 'JF017', 'JF019','JF020', 'AP083', 'AP084', 'JF021', 'JF022'}; % JF020 doesn't have histology yet. JF021 and JF022 to come, and many more
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
                                try
                                    AP_cellrasterJF({stimOn_times}, ...
                                        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                        trial_conditions(:, 4)});
                                catch
                                    try
                                        AP_cellrasterJF({stimOn_times}, ...
                                            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                            });
                                    catch
                                        AP_cellrasterJF({stimOn_times}, ...
                                            {stimIDs});
                                    end
                                end
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
