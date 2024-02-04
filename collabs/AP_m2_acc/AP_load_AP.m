
%% Find experiments with the task + cortical widefield + striatal ephys
close all;
myPaths;

animals={'AP100'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'AP_lcrGratingPassive'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

%% Load data from experiment 

curr_day = 1; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
date = thisDay;
verbose = false; % display load progress and some info figures
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true;

experiment = experiments(curr_day).experiment;
loadClusters = 0;

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
if ~isempty(qMetricsExist)

load(fullfile(savePath, 'qMetric.mat'))
load(fullfile(savePath, 'param.mat'))
bc_getQualityUnitType;
end
clearvars unitType 
AP_load_experiment;

curr_shank=NaN;
AP_cellraster({stimOn_times}, ...
    {stimIDs});

