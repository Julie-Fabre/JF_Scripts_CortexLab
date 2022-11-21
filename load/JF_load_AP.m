
%% Find experiments with the task + cortical widefield + striatal ephys
%close all;
%clear all;
myPaths;

animals={'JF070'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = ''; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

% allen_atlas_path = [allenAtlasPath, filesep, 'allenCCF/'];
% st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
% curr_plot_structure = find(strcmp(st.acronym, 'GPe'));
% histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
% load(histoFile)
% probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
% load(probe2ephysFile)
% min(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))
% max(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))

%% Load data from experiment 

curr_day = 5; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
date = thisDay;
verbose = false; % display load progress and some info figures
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true;

site = 1;%1,1; 2,4; 3,7
recording = []; 
% keep experiment with max n trials (= most likely not aborted error or end
% shank mapping) QQ change this in the future 
n_trials = zeros(size(experiments(curr_day).experiment,2), 1);
for iExperiment = experiments(curr_day).experiment
    [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiments(curr_day).experiment(iExperiment), 'block');
    load(block_filename)
    if isfield(block.events, 'stim_idValues')
        n_trials(iExperiment) = length(block.events.stim_idValues);
    elseif isfield(block.events, 'stimulusOnTimes')
        n_trials(iExperiment) = length(block.events.stimulusOnTimes);
    end
end 

experiment = experiments(curr_day).experiment(find(n_trials == max(n_trials)));
loadClusters = 0;
[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_includingCompressed',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = contains(ephysAPfile, '_g');
ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
% if ~isempty(qMetricsExist)
% 
% load(fullfile(savePath, 'qMetric.mat'))
% param = parquetread(fullfile(savePath, '_jf_parameters._jf_qMetrics.parquet'));
% bc_getQualityUnitType;
% end
clearvars unitType 
%load_parts.cam = true;
JF_load_experiment;
curr_shank=NaN;

trial_conditions(ismember(trial_conditions(:,1), [4]),1) = 4; % go 1
trial_conditions(ismember(trial_conditions(:,1), [7]),1) = 7; % go 2
trial_conditions(ismember(trial_conditions(:,1), [6]),1) = 10; % no go
trial_conditions(~ismember(trial_conditions(:,1), [4,7,10]),1) = 1;
%thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions, [4,7,10]) & trial_conditions(:,2)~=90;
thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions, [4,7,10]);



AP_cellrasterJF({stimOn_times(thisIndex), stimOn_times(thisIndex), stimOn_times(thisIndex)}, ...
    {trial_conditions(thisIndex,1), trial_conditions(thisIndex,2),...
    (trial_conditions(thisIndex,2)/-90)+(trial_conditions(thisIndex,1))});
% 
% AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))',stimOn_times}, ...
%     {trial_conditions(:,1),trial_conditions(:,2), ...
%     trial_conditions(:,3), movement_after200ms_and_type});
AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))'}, ...
    {trial_conditions(:,1),trial_conditions(:,2), ...
    trial_conditions(:,3)});




