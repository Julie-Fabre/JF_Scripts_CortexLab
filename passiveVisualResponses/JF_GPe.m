animal = 'AP084';
histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
load(histoFile)
probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
load(probe2ephysFile)
probe2ephys(5)=[];
recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
recInfo = readtable(recInfoFile);
allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);


protocolStrings = {'rating', 'ocation', 'nat'};
locations = {'GPe', 'Pallidum'};
%load probe
uniqueD = unique(recInfo.Date);
uniqueD = uniqueD(contains(uniqueD, '/')); %keep only dates (ie have a /) QQ
[dn,idx] = sort(datenum(uniqueD, 'dd/mm/yyyy'), 1, 'ascend');%sort by ascending order 
uniqueD = uniqueD(idx); 

ephysData = struct;
thisCount = 1;

for iProbe = 2:4 %probe

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
        a=recInfo.Site;
        b=cellfun(@numel,recInfo.Site);
        a(b>1 | b==0)={'0'};
        c=str2num(cell2mat(a));
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
            if numel(experiment > 1)
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
                thisCount = thisCount + 1;
            end
        end
    end
end