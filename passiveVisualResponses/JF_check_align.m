
%% load simple, check aligning
clear all;
animal = 'JF022';
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

ephysData = struct;
thisCount = 1;


for iProbe = 1:size(probe2ephys, 2)
    curr_day = probe2ephys(iProbe).day;
    curr_site = probe2ephys(iProbe).site;
    if isa(uniqueD, 'cell')
        cday = uniqueD{curr_day};
        dayBreaks = find(cday == '/' | cday == '-') ;
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
            | ~isnan(recInfo.timeline_not_started) ;
    elseif isa(recInfo.ephys_not_started, 'cell') && isa(recInfo.timeline_not_started, 'double') 
        theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
            | ~isnan(recInfo.timeline_not_started) ;
    elseif isa(recInfo.ephys_not_started, 'double') && isa(recInfo.timeline_not_started, 'string') 
        theseExpErrors = ~isnan(recInfo.ephys_not_started) ...
            | strlength(recInfo.timeline_not_started) >0 ;
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
    experiment=recInfo.Protocol_number(experiment); 
    verbose = false; % display load progress and some info figures
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;

    lfp_channel = 'all';
    close all;
    AP_load_experimentJF_WIP;
    if dontAnalyze==1%try another method
        AP_load_experimentJF;
    end
    % plot alignement stuffs

    % align to stim onset
    AP_cellrasterJF({stimOn_times}, {stimIDs})

    % align to movement onset
    AP_cellrasterJF({wheel_move_time}, {})
    
end