
%% Find experiments with the task
myPaths;
animal = 'JF086'; % mouse 
protocol = 'Grating'; % (this is the name of the Signals protocol) protcols live here: zserver/Code/Rigging/ExpDefinitions/JulieF/passiveWorld/visualResponsesNaiveBgg
%protocol = ''; % find all experiments, regardless of expdef name 
experiments = AP_find_experimentsJF(animal, protocol, true); % find experiments that contain protocol string 
experiments = experiments([experiments.ephys]); % QQ uncomment me  % only keep ephys experiments 

%% Load data from experiment 

curr_day = 2; % (set which day to use)
day = experiments(curr_day).day; % date
verbose = false; % display load progress and some info figures
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true; % QQ change me 

site = 1;%1,1; 2,4; 3,7
recording = []; 
experiment = experiments(curr_day).experiment(1);

loadClusters = 0; % whether to load keep 'good units' as defined by quality metrics or manual curation, if those exist
[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_includingCompressed',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = true; 
if isSpikeGlx
     [ephysKSfile,~] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys',site,recording);
    if isempty(dir([ephysKSfile filesep 'sync.mat']))
        warning('NO SYNC')
        %sync = syncFT(ephysAPfile, 385, ephysKSfile);
        if length(unique(sync)) == 1
            warning('check is sync is saved properly!')
        end
    end
end

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));

AP_load_experimentJF;
curr_shank=NaN;

AP_cellrasterJF({stimOn_times, stimOn_times, stimOn_times}, ...
    {trial_conditions(1:size(stimOn_times),2), trial_conditions(1:size(stimOn_times),3),...
    trial_conditions(1:size(stimOn_times),2) + trial_conditions(1:size(stimOn_times),3)});

