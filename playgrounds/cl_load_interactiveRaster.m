
%% Find experiments with the task + cortical widefield + striatal ephys
cl_myPaths;

animals = {'JF123'};
curr_animal = 1; % (set which animal to use)
animal = animals{curr_animal};

protocol = ''; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
experiments = experiments([experiments.ephys]);


%% Load data from experiment

curr_day = 8; % (set which day to use)
thisDate = experiments(curr_day).thisDate; % date


site = 1; %1,1; 2,4; 3,7
recording = [];

experiment = experiments(curr_day).experiment(2); %experiments(curr_day).experiment(1);%find(n_trials == max(n_trials));

[ephysAPfile, aa] = cl_cortexlab_filename(animal, thisDate, experiment, 'ephys_includingCompressed', site, recording);
if size(ephysAPfile, 2) == 2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end

ephysDirPath = cl_cortexlab_filename(animal, thisDate, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));


verbose = false; % display load progress and some info figures
load_parts.cam = false;%true;
load_parts.imaging = false;
load_parts.ephys = true;
loadClusters = 0;
cl_load_experiment;
curr_shank = NaN;
%AP_cellrasterJF({stimOn_times, stimOn_times}, {trial_conditions(:,1), trial_conditions(:,2)})


trial_conditions_clean = trial_conditions;
trial_conditions_clean(trial_conditions(:, 1) == 4, 1) = 100; %go1
trial_conditions_clean(trial_conditions(:, 1) == 12, 1) = 101; %go2
trial_conditions_clean(trial_conditions(:, 1) == 6, 1) = 102; %no go
%trial_conditions_clean(trial_conditions(:,1)==11,1) = 103;%no like
%trial_conditions_clean(trial_conditions(:,1)==13,1) = 104;%go like
theseImages = [100:102];


theseImages_trials = ismember(trial_conditions_clean(:, 1), theseImages) & ismember(trial_conditions_clean(:, 2), [-90]);
cl_cellraster({stimOn_times(theseImages_trials), stimOn_times(theseImages_trials), stimOn_times(theseImages_trials)}, ...
    {trial_conditions_clean(theseImages_trials, 1), ...
    trial_conditions_clean(theseImages_trials, 2), -89 + trial_conditions_clean(theseImages_trials, 1) + abs(trial_conditions_clean(theseImages_trials, 2))})

%% task
experiment = experiments(curr_day).experiment(1); %experiments(curr_day).experiment(1);%find(n_trials == max(n_trials));

[ephysAPfile, aa] = cl_cortexlab_filename(animal, thisDate, experiment, 'ephys_includingCompressed', site, recording);
if size(ephysAPfile, 2) == 2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end

ephysDirPath = cl_cortexlab_filename(animal, thisDate, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));


verbose = false; % display load progress and some info figures
load_parts.cam = true;
load_parts.imaging = false;
load_parts.ephys = true;
loadClusters = 0;
cl_load_experiment;
curr_shank = NaN;
%AP_cellrasterJF({stimOn_times, stimOn_times}, {trial_conditions(:,1), trial_conditions(:,2)})


AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(1:size(stimOn_times,1))'}, ...
     {trial_conditions(:,1),trial_conditions(:,2), ...
     trial_conditions(:,3)});