%% postdoc talk master 


%% 1. Naive 
% load data 
clear all;
loadVids = 0;
[nat_passive, session_data] =  cl_loadPerStimulusData('passive', 6, loadVids); 
save('/home/julie/Dropbox/MATLAB/nat_passive_full.mat', '-struct', 'nat_passive', '-v7.3');
close all;

clear all;
loadVids = 0;
[gratings_passive, session_data] = cl_loadPerStimulusData('passive', 1, loadVids);  
save('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat', '-struct', 'gratings_passive', '-v7.3');
close all;

clear all; % bug 
loadVids = 0;
locations_passive = cl_loadPerStimulusData('passive', 2, loadVids); 
save('/home/julie/Dropbox/MATLAB/locations_passive_full.mat', '-struct', 'locations_passive', '-v7.3');
close all;

% example cells

% population response 
passive_data = load('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat');
index = 1;
cl_population_PSTH;

% fraction cells 
cl_percentage_cells;


% selectivity: ex cells 

% selectivity: population 

% selectivity index 


%% 2. Task 


%% 3. Trained vs task 

%% All locations 
figure_num = 1;
passive_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat');
index = 5;
regions = {'CP', 'GPe', 'SNr'};
cl_location_PSTH;

figure_num = 2;
passive_data = load('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat');
index = 1;
regions = {'CP', 'GPe', 'SNr'};
cl_location_PSTH;

% location go go go 
figure_num = 3;
passive_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
index = 2;
regions = {'CP', 'GPe', 'SNr'};
cl_location_PSTH;

% location go no go
figure_num = 4;
passive_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
index = 2;
regions = {'CP', 'GPe', 'SNr'};
cl_location_PSTH;

% locations combined go go go and go no go 
figure_num = 5;
passive_data1 = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
passive_data2 = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
passive_data = cl_combineDataSets(passive_data1, passive_data2, 2, 2);
index = 1;
regions = {'CP', 'GPe', 'SNr'};
cl_location_PSTH;