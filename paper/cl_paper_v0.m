%% ~~ Naive ~~

%% Load data 
load_type = 'naive';
loadVids = 0;
for iExperimentType = 1:6
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/naive_data' num2str(iExperimentType) '.mat'], '-struct', 'task_data', '-v7.3');
end

load_type = 'taskGo';
loadVids = 0;
for iExperimentType = 1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/gogogo_data' num2str(iExperimentType) '.mat'], '-struct', 'task_data', '-v7.3');
end

load_type = 'TaskNoGo';
loadVids = 0;
for iExperimentType = 1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/gonogo_data' num2str(iExperimentType) '.mat'], '-struct', 'task_data', '-v7.3');
end

%% Fig 1: 
% Example cells visual 

% Percentage cells visual 

% Population cells visual 

% Example cells selectivity 

% Cell types 

%% Fig 2: selective cells, and throughout 

% Selectivity index (+ percentage "selective") 

% Population selectivity 

% Location 

%% ~~ Tasks - behavior data ~~

%% ~~ Tasks - neural data ~~

%% ~~ Suppl.: motion, movement ect. ; cell types ; delineating subregions ~~