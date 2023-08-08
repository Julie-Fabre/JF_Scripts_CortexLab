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
cl_plot_average_task; 


task_data = cl_loadPerStimulusData('task', 2); % 2 = 
gogogo = 0;
passive = 0;
cl_plot_average_task; 

gogogo = 0;
passive = 0;
cl_PC_analysis;

task_data_gogogo = cl_loadPerStimulusData('taskGo', 2); 
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
task_data = cl_loadPerStimulusData('task', 2); 
save('task_data_goNogo3.mat', '-struct', 'task_data');
task_data_gogogo = cl_loadPerStimulusData('taskGo', 2); 
save('task_data_gogogo.mat', '-struct', 'task_data_gogogo');
