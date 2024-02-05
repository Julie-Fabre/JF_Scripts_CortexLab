%% ~~ Naive ~~

%% Load data 
load_type = 'naive';
loadVids = 0;
for iExperimentType = 5:6% 1:6
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/naive_data' num2str(iExperimentType) '.mat'], '-struct', 'expData', '-v7.3');
end
% 5 (choiceworld!): Error using signrank
% No data remaining after removal of NaNs.
% 
% Error in cl_loadPerStimulusData (line 297)
%                         pvals_per_cond(iCond) = signrank(nanmean(curr_raster(logical(trial_cond_idx_single), 1:150), 2), ...
% 54/273

close all;
load_type = 'taskGo';
loadVids = 0;
for iExperimentType = 2%1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/gogogo_data' num2str(iExperimentType) '.mat'], '-struct', 'expData', '-v7.3');
end

close all;
load_type = 'taskNoGo';
loadVids = 0;
for iExperimentType = 2%1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/goNogo_data' num2str(iExperimentType) '.mat'], '-struct', 'expData', '-v7.3');
end


%% Fig 1: 
% Example cells visual 

% Percentage cells visual 

% Population cells visual 

% Example cells selectivity 

% Cell types 
% fast responses in GPe/SNr ? 
%% Fig 2: selective cells, and throughout 

% Selectivity index (+ percentage "selective") 

% Population selectivity 

% Location 

%% ~~ Tasks - behavior data ~~

%% ~~ Tasks - neural data ~~

%% ~~ Suppl.: motion, movement ect. ; cell types ; delineating subregions ~~

%% Suppl: histo
%% Suppl: bombcell 


%% others: RPE, ect. 