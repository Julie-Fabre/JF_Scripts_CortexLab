
%% Find experiments with the task
myPaths;
animal = 'JF086'; % mouse 
protocol = 'rating'; % (this is the name of the Signals protocol) protocols live here: zserver/Code/Rigging/ExpDefinitions/JulieF/passiveWorld/visualResponsesNaiveBgg
%protocol = ''; % find all experiments, regardless of expdef name 
experiments = AP_find_experimentsJF(animal, protocol, true); % find experiments that contain protocol string 
experiments = experiments([experiments.ephys]); % only keep ephys experiments 

%% Load data from experiment 

curr_day = 1; % (set which day to use)
day = experiments(curr_day).day; % date
verbose = false; % display load progress and some info figures
load_parts.cam=false; % load eye and face videos, align to timeline time
load_parts.imaging=false; % this is if we had widefield data
load_parts.ephys=true; % load kilosorted data

site = 1;
recording = []; 
experiment = experiments(curr_day).experiment(site);

loadClusters = 0; % whether to load keep 'good units' as defined by quality metrics or manual curation, if those exist
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_includingCompressed',site,recording);
if size(ephysAPfile,2) == 2 %keep only AP file, not any meta or lfp
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = true; % always true in your case
ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));

JF_load_experiment;
curr_shank=NaN;

%% look at data
% This is a raster/PSTH viewer to explore cells, for a given unit it'll
% plot all spikes aligned to different times (first input) and can split
% those into different conditions (second input). The lefthand plot shows
% one dot for each unit (depth by rate) which are clickable, the center
% plot shows the template waveform , the righthand plot shows the spikes,
% the bottom plot shows the template amplitude over time (to see if it
% drifts, etc)
%
% controls:
% * left/right arrows switch between stim, move, and outcome aligned
% * pageup/pagedown combines all trials or splits stim, move directions,
% and reward/punish, respectively
% * you can press 'm' and define borders by clicking on the lefthand plot
% to group cells into multiunit

AP_cellrasterJF({stimOn_times, stimOn_times, stimOn_times}, ...
    {trial_conditions(1:size(stimOn_times),2), trial_conditions(1:size(stimOn_times),3),...
    trial_conditions(1:size(stimOn_times),2) + trial_conditions(1:size(stimOn_times),3)});

