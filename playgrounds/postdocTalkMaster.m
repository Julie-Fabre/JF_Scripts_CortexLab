%% postdoc talk master 

%% ~~ 0. load all data ~~
loadVids = 0;
[task_data, session_data] = cl_loadPerStimulusData('task', 2, loadVids); 
%save('/home/julie/Dropbox/MATLAB/task_data_goNogo_session.mat', '-struct', 'session_data', '-v7.3');
save('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat', '-struct', 'task_data', '-v7.3');
close all;

clear all;
loadVids = 0;
[task_data_gogogo, session_data] = cl_loadPerStimulusData('taskGo', 2, loadVids); 
%save('/home/julie/Dropbox/MATLAB/task_data_gogogo_session.mat', '-struct', 'session_data', '-v7.3');
save('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat', '-struct', 'task_data_gogogo', '-v7.3');
close all;

clear all;
loadVids = 0;
[task_data_passive, session_data] = cl_loadPerStimulusData('passive', 5, loadVids); 
%save('/home/julie/Dropbox/MATLAB/task_data_passive_session.mat', '-struct', 'session_data', '-v7.3');
save('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat', '-struct', 'task_data_passive', '-v7.3');
close all;

clear all;
loadVids = 0;
[nat_passive, session_data] =  cl_loadPerStimulusData('passive', 6, loadVids); 
save('/home/julie/Dropbox/MATLAB/nat_passive_full.mat', '-struct', 'nat_passive', '-v7.3');
close all;

clear all;
loadVids = 0;
[gratings_passive, session_data] = cl_loadPerStimulusData('passive', 1, loadVids);  
save('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat', '-struct', 'gratings_passive', '-v7.3');

clear all;
loadVids = 0;
locations_passive = cl_loadPerStimulusData('passive', 2, loadVids); 
save('/home/julie/Dropbox/MATLAB/locations_passive_full.mat', '-struct', 'locations_passive', '-v7.3');

clear all;
loadVids = 0;
task_data_passive2 = cl_loadPerStimulusData('passive', 4, loadVids); 
save('/home/julie/Dropbox/MATLAB/task_data_passive2.mat', '-struct', 'task_data_passive2', '-v7.3');

clear all;
loadVids = 0;
[nat_passive, session_data] = cl_loadPerStimulusData('passive', 3, loadVids); 
%save('/home/julie/Dropbox/MATLAB/task_data_passive_session.mat', '-struct', 'session_data', '-v7.3');
save('/home/julie/Dropbox/MATLAB/nat_passive_full2.mat', '-struct', 'nat_passive', '-v7.3');
close all;


%% ~~ 1. Naive ~~

%% example cells

%% cell types 
% waveform duration, post spike suppression, prop ISI, firing rate 
passive_data1 = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
passive_data2 = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
passive_data = cl_combineDataSets(passive_data1, passive_data2, 2, 2);
passive_data3 = load('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat');
passive_data = cl_combineDataSets(passive_data, passive_data3, 1, 1);
cl_cellTypes;

%% population response 
passive_data = load('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat');
index = 1;
cl_population_PSTH;


%% fraction cells 
cl_percentage_cells;


%% selectivity: ex cells 

%% selectivity: population 

%% selectivity index 

regions = {'CP', 'GPe', 'SNr'};
cl_selectivity;

%% ~~ 2. Task ~~


%% ~~ 3. Trained vs task ~~
close all;
cl_plot_PSTHs(1, 0);%contra
cl_plot_PSTHs(0, 1);%center


subplot(3,3,1)
clim([-5,5])
subplot(3,3,2)
clim([-5,5])
subplot(3,3,3)
clim([-5,5])

subplot(3,3,4)
clim([-1,1])
subplot(3,3,5)
clim([-1,1])
subplot(3,3,6)
clim([-1,1])

subplot(3,3,7)
clim([-1.5,1.5])
subplot(3,3,8)
clim([-1.5,1.5])
subplot(3,3,9)
clim([-1.5,1.5])
prettify_plot;
%% Trained: perc. responsive
cl_fraction_cells_summary;
%cl_visResp_summarize;


%% Trained: dprime 
cl_dprime_summarize;
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