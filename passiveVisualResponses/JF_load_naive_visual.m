
%% load naive visual recordings
myPaths;
corona = 0;
mice = {'AP080', 'AP082', 'AP083', 'AP084', 'JF017', 'JF019', 'JF020', 'JF022', 'JF023',...
    'JF024', 'JF030', 'JF045','JF047','JF048','JF049','JF052'}; % JF020 doesn't have histology yet. JF021 and JF022 to come, and many more
locations = {'CP', 'GPe', 'GPi', 'STN','SNr'};
protocols = {'JF_natural_images', 'JF_locations', 'JF_GratingsPassive'};
ephysDataEnd = struct;
recordingInfo = readtable('/home/julie/Dropbox/Analysis/RecNew.csv');
allenAt = loadStructureTreeJF([allenAtlasPath filesep 'allenCCF/structure_tree_safe_2017.csv']);
allen_atlas_path = [allenAtlasPath filesep 'allenCCF/'];
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);

thisCount = 1;
for iMouse = 13:size(mice, 2)
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
            close all;
            corrSite = find(arrayfun(@(x) probe2ephys(x).site == iSite, 1:numel(probe2ephys)));
            probe = intersect(corrDay, corrSite);
            if length(probe)>2
                allprobe = probe;
                %probe = allprobe(iProbe);
                for iProbe = 1:size(allprobe,2)
                    probe = allprobe(iProbe);
                    thisP = probe;
                    for iLocation = 1:size(locations, 2)


                    %disp(probe)
                    %thisP = probe;
                    if ~isempty(probe)
                        this_ccf = probe_ccf(probe); %find if contains location of interest
                        if ~any(structfun(@isempty, this_ccf)) %check hewre
                            curr_plot_structure = find(strcmp(st.acronym, locations{iLocation}));
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
                                if isfield(recInfo, 'error')
                                    if iscell(recInfo.error(1))
                                        errorOK = cellfun(@isempty, recInfo.error(thisProtocol));
                                    else
                                        errorOK = ones(length(experimentThese), 1);
                                    end
                                else
                                    errorOK = ones(length(experimentThese), 1);
                                end
                                experimentThese = experimentThese(ephysOK & timelineOK & errorOK);


                                %experiment = experiment(probe2ephys(probe).site);
                                verbose = false; % display load progress and some info figures
                                load_parts.cam = false;
                                load_parts.imaging = false;
                                load_parts.ephys = true;
                                site = probe2ephys(probe).site;
                                for iExperiment = 1:size(experimentThese, 2)
                                    experiment = experimentThese(iExperiment);
                                    loadClusters = 0;

                                    isSpikeGlx = 0;
                                    if sum(animal == 'JF045')==5 && iSite == 3%hacky 
                                        recording = 1;
                                    else
                                        recording = [];
                                    end
                                    loadLFP=0;
                                    try
                                            recording = [];
                                        [ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap',site,recording);
                                       isSpikeGlx=1;
                                    AP_load_experimentJF;
                                    close all;
                                    if isfield(probe2ephys, 'shank')
                                        if ~isnan(probe2ephys(thisP).shank)
                                            curr_shank = probe2ephys(thisP).shank;
                                            theseChannelPositions = [(curr_shank-1) * 250, (curr_shank-1)*250 + 32];
                                            theseChannels = ismember(channel_positions(:,1), theseChannelPositions);
                                            theseTemplates = ismember(template_xdepths, theseChannelPositions);
                                            theseSpikes = ismember(spike_xdepths, theseChannelPositions);
                                            spike_times = spike_times(theseSpikes);
                                            spike_templates = spike_templates(theseSpikes);
                                            %rename 
                                            good_templates_idx = unique(spike_templates);
                                            new_spike_idx = nan(max(spike_templates), 1);
                                            new_spike_idx(good_templates_idx) = 1:length(good_templates_idx);
                                            spike_templates = new_spike_idx(spike_templates);

                                            template_amplitudes = template_amplitudes(theseSpikes);
                                            template_depths = template_depths(theseTemplates);
                                            channel_positions = channel_positions(theseChannels,:);
                                            templates = templates(theseTemplates,:,theseChannels); 
                                        end
                                    end
                                    %                                 try
                                    %                                     AP_cellrasterJF({stimOn_times}, ...
                                    %                                         {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                    %                                         trial_conditions(:, 4)});
                                    %                                 catch
                                    %                                     try
                                    %                                         AP_cellrasterJF({stimOn_times}, ...
                                    %                                             {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                    %                                             });
                                    %                                     catch
                                    %                                         AP_cellrasterJF({stimOn_times}, ...
                                    %                                             {stimIDs});
                                    %                                     end
                                    %                                 end
                                    %template_depths = template_depths; 
                                    theseTemplates = find(template_depths >= min(theseDepths) & template_depths <= max(theseDepths)); %correct depth units
                                    %find closest depth and get closest location
    %                                 [~,dv_sort_idx] = sort(this_ccf.trajectory_coords(:,2));
    % 
    %                                 probe_trajectory_depths = ...
    %                                     pdist2(this_ccf.trajectory_coords, ...
    %                                     this_ccf.trajectory_coords((dv_sort_idx == 1),:))*10;
                                    thisStructOnly = find(this_ccf.trajectory_areas==curr_plot_structure);
                                    A = repmat(this_ccf.probe_depths, [1, length(template_depths(theseTemplates))]);
                                    [minValue, closestIndex] = min(abs(A-template_depths(theseTemplates)'));
                                    closestLocation = this_ccf.trajectory_coords(closestIndex, :);
                                    thisRegion = this_ccf.trajectory_areas(closestIndex);
                                    %raster general and for each trial condition
                                    theseSpikes = ismember(spike_templates, theseTemplates);
                                    %save info
                                    %if bad_flipper == 0
                                    ephysDataEnd(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
                                    ephysDataEnd(thisCount).spike_templates = spike_templates(theseSpikes);
                                    try
                                        ephysDataEnd(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
                                    end
                                    ephysDataEnd(thisCount).templatesOI = theseTemplates;
                                    ephysDataEnd(thisCount).template_depths = template_depths(theseTemplates);

                                    ephysDataEnd(thisCount).template_wf = waveforms(theseTemplates, :);

                                    ephysDataEnd(thisCount).template_location = closestLocation;
                                    ephysDataEnd(thisCount).protocol = expDef;
                                    ephysDataEnd(thisCount).animal = animal;
                                    ephysDataEnd(thisCount).date = day;
                                    ephysDataEnd(thisCount).curr_day = curr_day;
                                    ephysDataEnd(thisCount).probe = thisP;
                                    ephysDataEnd(thisCount).site = site;
                                    ephysDataEnd(thisCount).experiment = experiment;
                                    ephysDataEnd(thisCount).location = locations(iLocation);
                                    ephysDataEnd(thisCount).stimOn_times = stimOn_times;
                                    ephysDataEnd(thisCount).stimIDs = stimIDs;
                                    if isfield(probe2ephys, 'shank')
                                        ephysDataEnd(thisCount).shank = probe2ephys(thisP).shank;
                                    end
                                   % if bad_flipper == 0
                                    %ephysDataEnd(thisCount).template_label = template_label(theseTemplates);
                                    %end
                                    %ephysDataEnd(thisCount).template_label_good = template_label_good(theseTemplates);
                                    thisCount = thisCount + 1;
                                    catch
                                    end
                                end
                            end
                        end
                    end
                    end
                end
            else
            for iLocation = 1:size(locations, 2)
                
                
                %disp(probe)
                thisP = probe;
                if ~isempty(probe)
                    this_ccf = probe_ccf(probe); %find if contains location of interest
                    if ~any(structfun(@isempty, this_ccf)) %check hewre
                        curr_plot_structure = find(strcmp(st.acronym, locations{iLocation}));
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
                            if isfield(recInfo, 'error')
                                if iscell(recInfo.error(1))
                                    errorOK = cellfun(@isempty, recInfo.error(thisProtocol));
                                else
                                    errorOK = ones(length(experimentThese), 1);
                                end
                            else
                                errorOK = ones(length(experimentThese), 1);
                            end
                            experimentThese = experimentThese(ephysOK & timelineOK & errorOK);


                            %experiment = experiment(probe2ephys(probe).site);
                            verbose = false; % display load progress and some info figures
                            load_parts.cam = false;
                            load_parts.imaging = false;
                            load_parts.ephys = true;
                            site = probe2ephys(probe).site;
                            for iExperiment = 1:size(experimentThese, 2)
                                experiment = experimentThese(iExperiment);
                                loadClusters = 0;

                                isSpikeGlx = 0;
                                if sum(animal == 'JF045')==5 && iSite == 3%hacky 
                                    recording = 1;
                                else
                                    recording = [];
                                end
                                loadLFP=0;
                                try
                                AP_load_experimentJF;
                                close all;
                                if isfield(probe2ephys, 'shank')
                                    if ~isnan(probe2ephys(thisP).shank)
                                        theseChannelPositions = [(curr_shank-1) * 250, (curr_shank-1)*250 + 32];
                                        theseChannels = ismember(channel_positions(:,1), theseChannelPositions);
                                        theseTemplates = ismember(template_xdepths, theseChannelPositions);
                                        theseSpikes = ismember(spike_xdepths, theseChannelPositions);
                                        spike_times = spike_times(theseSpikes);
                                        spike_templates = spike_templates(theseSpikes);
                                        %rename 
                                        good_templates_idx = unique(spike_templates);
                                        new_spike_idx = nan(max(spike_templates), 1);
                                        new_spike_idx(good_templates_idx) = 1:length(good_templates_idx);
                                        spike_templates = new_spike_idx(spike_templates);

                                        template_amplitudes = template_amplitudes(theseSpikes);
                                        template_depths = template_depths(theseTemplates);
                                        channel_positions = channel_positions(theseChannels,:);
                                        templates = templates(theseTemplates,:,theseChannels); 
                                    end
                                end
                                %                                 try
                                %                                     AP_cellrasterJF({stimOn_times}, ...
                                %                                         {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                %                                         trial_conditions(:, 4)});
                                %                                 catch
                                %                                     try
                                %                                         AP_cellrasterJF({stimOn_times}, ...
                                %                                             {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
                                %                                             });
                                %                                     catch
                                %                                         AP_cellrasterJF({stimOn_times}, ...
                                %                                             {stimIDs});
                                %                                     end
                                %                                 end
                                
                                theseTemplates = find(template_depths >= min(theseDepths) & template_depths <= max(theseDepths)); %correct depth units
                                %find closest depth and get closest location
%                                 [~,dv_sort_idx] = sort(this_ccf.trajectory_coords(:,2));
% 
%                                 probe_trajectory_depths = ...
%                                     pdist2(this_ccf.trajectory_coords, ...
%                                     this_ccf.trajectory_coords((dv_sort_idx == 1),:))*10;
                                thisStructOnly = find(this_ccf.trajectory_areas==curr_plot_structure);
                                A = repmat(this_ccf.probe_depths, [1, length(template_depths(theseTemplates))]);
                                [minValue, closestIndex] = min(abs(A-template_depths(theseTemplates)'));
                                closestLocation = this_ccf.trajectory_coords(closestIndex, :);
                                thisRegion = this_ccf.trajectory_areas(closestIndex);
                                %raster general and for each trial condition
                                theseSpikes = ismember(spike_templates, theseTemplates);
                                %save info
                                %if bad_flipper == 0
                                ephysDataEnd(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
                                ephysDataEnd(thisCount).spike_templates = spike_templates(theseSpikes);
                                try
                                    ephysDataEnd(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
                                end
                                ephysDataEnd(thisCount).templatesOI = theseTemplates;
                                ephysDataEnd(thisCount).template_depths = template_depths(theseTemplates);

                                ephysDataEnd(thisCount).template_wf = waveforms(theseTemplates, :);

                                ephysDataEnd(thisCount).template_location = closestLocation;
                                ephysDataEnd(thisCount).protocol = expDef;
                                ephysDataEnd(thisCount).animal = animal;
                                ephysDataEnd(thisCount).date = day;
                                ephysDataEnd(thisCount).curr_day = curr_day;
                                ephysDataEnd(thisCount).probe = thisP;
                                ephysDataEnd(thisCount).site = site;
                                ephysDataEnd(thisCount).experiment = experiment;
                                ephysDataEnd(thisCount).location = locations(iLocation);
                                ephysDataEnd(thisCount).stimOn_times = stimOn_times;
                                ephysDataEnd(thisCount).stimIDs = stimIDs;
                                if isfield(probe2ephys, 'bank')
                                    ephysDataEnd(thisCount).bank = probe2ephys(thisP).bank;
                                end
                               % if bad_flipper == 0
                                %ephysDataEnd(thisCount).template_label = template_label(theseTemplates);
                                %end
                                %ephysDataEnd(thisCount).template_label_good = template_label_good(theseTemplates);
                                thisCount = thisCount + 1;
                                catch
                                end
                            end
                        end
                    end
                end
            end
            end
        end
    end
end


