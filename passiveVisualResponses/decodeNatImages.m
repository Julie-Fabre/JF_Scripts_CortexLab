
%% Natural images decoding

%% 1. load all data (DMS, PS), check alignement (plot, + look at alignement to movement) -> dont manually
toLoad = readtable('C:/Users/Julie/Dropbox/NaturalImagesDatasetsForDecodingGrant.ods'); % get which animals and probes to load from each animal, and depth corr if necessary

theseAnimals = unique(toLoad.Animal);
ephysData = struct;
thisCount = 1;
for iAnimal = 1:size(theseAnimals)

    animal = theseAnimals{iAnimal};
    theseProbes = toLoad.Probe(find(strcmp(toLoad.Animal, animal))); %get probes to load
    theseDepthsCorr = toLoad.Depth(find(strcmp(toLoad.Animal, animal))); %get probes to load
    theseLocationsCorr = toLoad.Location(find(strcmp(toLoad.Animal, animal))); %get probes to load

    protString = 'mage';
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)
    recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
    recInfo = readtable(recInfoFile);
    allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);

    locations = {'CP'};
    %load probe
    uniqueD = unique(recInfo.Date);
    if isa(uniqueD, 'datetime')
        [dn, idx] = sort(uniqueD, 1, 'ascend');
    else
        uniqueD = uniqueD(contains(uniqueD, '/') | contains(uniqueD, '-')); %keep only dates (ie have a /) QQ
        try
            [dn, idx] = sort(datenum(uniqueD, 'dd/mm/yyyy'), 1, 'ascend'); %sort by ascending order
        catch
            [dn, idx] = sort(datenum(uniqueD, 'dd-mm-yyyy'), 1, 'ascend'); %sort by ascending order
        end
    end
    uniqueD = uniqueD(idx);

    for iP = 1:length(theseProbes)
        iProbe = theseProbes(iP);
        curr_day = probe2ephys(iProbe).day;
        curr_site = probe2ephys(iProbe).site;
        if isa(uniqueD, 'cell')
            cday = uniqueD{curr_day};
            dayBreaks = find(cday == '/' | cday == '-');
            if dayBreaks(1) == 2 && dayBreaks(2) == 4
                cday = strcat(cday(dayBreaks(2)+1:end), '-0', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-0', cday(1:dayBreaks(1)-1)); %in standard format
            elseif dayBreaks(1) == 2
                cday = strcat(cday(dayBreaks(2)+1:end), '-', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-0', cday(1:dayBreaks(1)-1)); %in standard format
            elseif dayBreaks(2) == 5
                cday = strcat(cday(dayBreaks(2)+1:end), '-0', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-', cday(1:dayBreaks(1)-1)); %in standard format
            else
                cday = strcat(cday(dayBreaks(2)+1:end), '-', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-', cday(1:dayBreaks(1)-1)); %in standard format
            end
            theseExpDates = strcmp(recInfo.Date, uniqueD{curr_day});
        else
            cday = uniqueD(curr_day);
            day = cday;
            theseExpDates = isbetween(recInfo.Date, uniqueD(curr_day), uniqueD(curr_day));
        end


        if ~iscell(recInfo.Site)
            theseExpSites = recInfo.Site == curr_site;
        else
            a = recInfo.Site;
            b = cellfun(@numel, recInfo.Site);
            a(b > 1 | b == 0) = {'0'};
            c = str2num(cell2mat(a));
            theseExpSites = c == curr_site;
        end
        if exist('recInfo.error', 'var')
            theseExpErrors = contains(recInfo.error, 'es') | ~isnan(recInfo.ephys_not_started) ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'string') && isa(recInfo.timeline_not_started, 'string')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | strlength(recInfo.timeline_not_started) > 0;
        elseif isa(recInfo.ephys_not_started, 'double') && isa(recInfo.timeline_not_started, 'double')
            theseExpErrors = ~isnan(recInfo.ephys_not_started) ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'cell') && isa(recInfo.timeline_not_started, 'double')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'cell') && isa(recInfo.timeline_not_started, 'cell')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | strlength(recInfo.timeline_not_started) > 0;
        elseif isa(recInfo.ephys_not_started, 'double') && isa(recInfo.timeline_not_started, 'string')
            theseExpErrors = ~isnan(recInfo.ephys_not_started) ...
                | strlength(recInfo.timeline_not_started) > 0;
        end
        theseExp = find(theseExpDates & theseExpSites & ~theseExpErrors);
        theseProtocols = recInfo.Protocol(theseExp);

        experiment = theseExp(contains(theseProtocols, protString));
        protocol = 'atural';
        experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
        %experiments = experiments([experiments.imaging] & [experiments.ephys]);
        experiments = experiments([experiments.ephys]);

        site = recInfo.Site(experiment);

        day = experiments(curr_day).day; % date
        experiment = recInfo.Protocol_number(experiment);
        verbose = false; % display load progress and some info figures
        load_parts.cam = false;
        load_parts.imaging = false;
        load_parts.ephys = true;

        lfp_channel = 'all';
        close all;
        AP_load_experimentJF_WIP;
        if dontAnalyze == 1 %try another method
            AP_load_experimentJF;
        end
        % plot alignement stuffs
        if isempty(theseDepthsCorr{iP, :})
            this_ccf = probe_ccf(iProbe);
            theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
            theseLocationsInterest = contains(theseLocations, locations{1});
            theseDepths = this_ccf.probe_depths(theseLocationsInterest);
        else
            theseDepths_op = strsplit(theseDepthsCorr{iP, :}, '-');
            clearvars theseDepths
            theseDepths(1) = str2num(theseDepths_op{1, 1});
            theseDepths(2) = str2num(theseDepths_op{1, 2});
        end
        if ~isempty(theseDepths)
            %save data in structure
            theseTemplates = find(template_depths >= min(theseDepths) & template_depths <= max(theseDepths)); %correct depth units
            %find closest depth and get closest location
            if isempty(theseDepthsCorr{iP, :})
                A = repmat(this_ccf.probe_depths, [1, length(template_depths)]);
                [minValue, closestIndex] = min(abs(A-template_depths'));
                closestLocation = this_ccf.trajectory_coords(closestIndex, :);
            else
                closestLocation = NaN;
            end
            %raster general and for each trial condition
            theseSpikes = ismember(spike_templates, theseTemplates);
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
            ephysData(thisCount).location = theseLocationsCorr(iP);
            ephysData(thisCount).stimOn_times = stimOn_times;
            ephysData(thisCount).stimIDs = stimIDs;
            ephysData(thisCount).trial_conditions = trial_conditions;
            thisCount = thisCount + 1;
        end
        keep toLoad  theseAnimals  ephysData  thisCount iP iAnimal animal theseProbes theseDepthsCorr theseLocationsCorr protString ...
            probe_ccf probe2ephys recInfo allenAt locations uniqueD

    end
end
%save data on dropbox

%% 2. for each cell FR 1/2 trials against 1/2 -> r2, plot distribution

%% 3. one-way ANOVA: -way anova for each cell: is there a main effect "stim" on firing rate  ? pvalue for each cell. can then plot the districbution of pvalues and see which ones are signif.

%% 4. decoder
