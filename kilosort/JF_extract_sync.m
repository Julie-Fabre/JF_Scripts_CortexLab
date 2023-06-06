
%% JF_extract_sync(mouse)
%function JF_extract_sync(animal)
animals = {'JF093'};
for iAnimal = 1:length(animals)
animal = animals{iAnimal};
myPaths;
protocol = []; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
for iExperiment = 1:size(experiments, 1)
    %experiment = experiments(iExperiment);
    date = experiments(iExperiment).day;
    [ephysParentDir, ~] = AP_cortexlab_filenameJF(animal, date, [], 'ephys_parentDir', [], []);
    siteDir = dir([ephysParentDir, filesep, 'site*']);
    siteDir = siteDir(~contains({siteDir.name}, {'shank'}));
    for iSite = 1:size(siteDir, 1)
         site = str2num(siteDir(iSite).name(end));
        [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, date, [], 'ephys', site, []);
        if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
            [ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, [], 'ephys_includingCompressed', site, []);
            if contains(ephysAPfile, '.cbin')
                onlySaveSyncChannel= 1;
                decompDataFile = bc_extractCbinData(ephysAPfile, [], [], 0, ephysKSfile, onlySaveSyncChannel);

            else
            if size(ephysAPfile, 2) == 2 && iscell(ephysAPfile) %keep only ap
                ephysAPfile = ephysAPfile{1};
            end
            isSpikeGlx = contains(ephysAPfile, '_g'); %spike glx (2.0 probes) or open ephys (3A probes)?
            if isSpikeGlx
                if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
                    syncFT(ephysAPfile, 385, ephysKSfile);
                end
            else
                %AP_preprocess_phase3_newOEJF_onlysync(animal, date, site)
            end
            end
        
        end

    end
end
end
%end
