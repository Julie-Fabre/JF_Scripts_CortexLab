% animals, curr_animal, protocol, curr_day , site, experiment, recording,
% loadClisters load_parts.cam=false; load_parts.imaging=false; load_parts.ephys=true;

myPaths;

animals={'JF070'};
curr_animal = 1; % (set which animal to use)
animal = animals{curr_animal};

experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
date = thisDay;
verbose = false; % display load progress and some info figures


[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_ap',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = contains(ephysAPfile, '_g');
if isSpikeGlx
     [ephysKSfile,~] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys',site,recording);
    if isempty(dir([ephysKSfile filesep 'sync.mat']))
        warning('NO SYNC')
        %sync = syncFT(ephysAPfile, 385, ephysKSfile);
        if find(unique(sync)==0) < 10
            warning('check is sync is saved properly!')
        end
    end
end

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
if ~isempty(qMetricsExist)

    load(fullfile(savePath, 'qMetric.mat'))
    try
        param = parquetread(fullfile(savePath, '_jf_parameters._jf_qMetrics.parquet'));
        bc_getQualityUnitType;
        disp('loading quality metrics...')
    catch
        warning([animal, ' ', date ' - no quality metrics'])
    end
else
    warning([animal, ' ', date ' - no quality metrics'])
end
%clearvars unitType 

AP_load_experimentJF;





