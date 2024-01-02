%% load data
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
[task_data_passive, session_data] = cl_loadPerStimulusData('passive', 5, loadVids); 
%save('/home/julie/Dropbox/MATLAB/task_data_passive_session.mat', '-struct', 'session_data', '-v7.3');
save('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat', '-struct', 'task_data_passive', '-v7.3');
close all;

clear all;
[nat_passive, session_data] =  cl_loadPerStimulusData('passive', 6); 
%save('/home/julie/Dropbox/MATLAB/nat_passive_full.mat', '-struct', 'nat_passive', '-v7.3');
close all;

clear all;
[gratings_passive, session_data] = cl_loadPerStimulusData('passive', 1); 
%save('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat', '-struct', 'gratings_passive', '-v7.3');

clear all; % bug 
locations_passive = cl_loadPerStimulusData('passive', 2); 
save('/home/julie/Dropbox/MATLAB/locations_passive_full.mat', '-struct', 'locations_passive', '-v7.3');

clear all;
task_data_passive2 = cl_loadPerStimulusData('passive', 4); 
save('/home/julie/Dropbox/MATLAB/task_data_passive2.mat', '-struct', 'task_data_passive2', '-v7.3');

%% Naive: example cells 

%% Naive: population 

%% Naive: perc. responsive

%% Naive: selectivity

%% Naive: which stimuli respond? 

%% Naive: locations 

%% Trained: example cells, population
cl_average_response_per_image;

%% Trained: perc. responsive

%% Trained: dprime 
