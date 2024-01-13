%% passive_master script 

task_data = cl_loadAverageData('task'); 

% % do cells encode visual stimuli? 
% cl_loadAverageData; % check we don't have any duplicates 
% cl_population_PSTH;
% cl_location_PSTH;
% cl_percentage_cells;
% 
% % cell classification 
% cl_population_PSTH_celltype;
% cl_celltype_example;
% cl_GPe_celltype_playground; 
% 
% % cell selctivity/encoding of specific visual stimuli 
passive_data = cl_loadPerStimulusData('passive');
gogogo = 0;
passive = 1;
cl_plot_average_task; 


task_data = cl_loadPerStimulusData('task', 2); % 2 = 
gogogo = 0;
passive = 0;
cl_plot_average_task; 

gogogo = 0;
passive = 0;
cl_PC_analysis;

task_data_gogogo = cl_loadPerStimulusData('taskGo', 1); 

task_data_here = load('task_data_gogogo.mat');
gogogo = 1;
passive =0;
cl_plot_average_task; 

gogogo = 1;
passive =0;
cl_PC_analysis;
% cl_selectivity;
% %decoding? 
% cl_location_selectivity; 
% 
% % movement stuffs
% cl_motionIndex; 
% 
% % rastermap/umap


% DMS vs PS 
%% 
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

clear all; % bug 
loadVids = 0;
locations_passive = cl_loadPerStimulusData('passive', 2, loadVids); 
save('/home/julie/Dropbox/MATLAB/locations_passive_full.mat', '-struct', 'locations_passive', '-v7.3');


clear all;
loadVids = 0;
task_data_passive2 = cl_loadPerStimulusData('passive', 4, loadVids); 
save('/home/julie/Dropbox/MATLAB/task_data_passive2.mat', '-struct', 'task_data_passive2', '-v7.3');






%% D-prime 
cl_dprime_summary;

%% Average response to each image 
cl_average_response_per_image;

%% Fraction of cells 
cl_fraction_cells_summary;
cl_fraction_cells;