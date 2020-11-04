
%% load naive visual recordings
corona = 0;
mice = {'AP080', 'AP082'}; %, 'JF017'}; % JF017 doesn't have histology yet
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
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)

    for iDate = 1:size(experiments, 1)
        corrDay = find(arrayfun(@(x) probe2ephys(x).day == 2, 1:numel(probe2ephys)));
        for iCorrDay = 1:size(corrDay, 2)
            sites(iCorrDay) = probe2ephys(corrDay(iCorrDay)).site;
        end
        for iSite = 1:size(sites, 2)
            corrSite = find(arrayfun(@(x) probe2ephys(x).site == iSite, 1:numel(probe2ephys)));
            for iLocation = 1:size(locations, 2)
                probe = intersect(corrDay, corrSite);
                this_ccf = probe_ccf(probe); %find if contains location of interest
                if ~any( structfun(@isempty, this_ccf) )
                    theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
                    theseLocationsInterest = contains(theseLocations, locations{iLocation});
                    theseDepths = this_ccf.probe_depths(theseLocationsInterest);
                    if ~isempty(find(theseDepths))
                        %load experiment
                        curr_day = probe2ephys(probe).day;
                        day = experiments(curr_day).day;
                        experimentThese = experiments(curr_day).experiment; % experiment number
                        %experiment = experiment(probe2ephys(probe).site);
                        verbose = false; % display load progress and some info figures
                        load_parts.cam = false;
                        load_parts.imaging = false;
                        load_parts.ephys = true;
                        site = probe2ephys(probe).site;
                        for iExperiment = 1:size(experimentThese, 2)
                            experiment = experimentThese(iExperiment);
                            AP_load_experimentJF;
                            theseTemplates = template_depths >= min(theseDepths) & template_depths <= max(theseDepths); %correct depth units
                            %find closest depth and get closest location
                            A = repmat(this_ccf.probe_depths, [1, length(template_depths)]);
                            [minValue, closestIndex] = min(abs(A-template_depths'));
                            closestLocation = this_ccf.trajectory_coords(closestIndex, :);
                            %raster general and for each trial condition
                            theseSpikes = ismember(spike_templates, theseTemplates);
                            %save info
                            ephysData(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
                            ephysData(thisCount).spike_templates = spike_templates(theseSpikes);
                            ephysData(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
                            ephysData(thisCount).templatesOI = theseTemplates;
                            ephysData(thisCount).template_depths = template_depths(theseTemplates);
                            ephysData(thisCount).template_wf = waveforms(theseTemplates, waveform_peak(~isnan(new_spike_idx))- ...
                                4:waveform_peak(~isnan(new_spike_idx))+4);
                            ephysData(thisCount).template_location = closestLocation;
                            ephysData(thisCount).protocol = expDef;
                            ephysData(thisCount).animal = animal;
                            ephysData(thisCount).date = day;
                            ephysData(thisCount).site = site;
                            ephysData(thisCount).experiment = experiment;
                            ephysData(thisCount).location = locations(iLocation);
                            thisCount = thisCount + 1;
                        end
                    end
                end
            end
        end
    end
end
%     indexes = probe_ccf .* 10;
%     av(indexes)
%     av(round(probe_ccf(:, 1)), probe_ccf(:, 2), probe_ccf(:, 3))

%% plot responses