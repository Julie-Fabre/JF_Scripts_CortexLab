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
task_data = cl_loadPerStimulusData('task'); 
cl_plot_average_task; 
% cl_selectivity;
% %decoding? 
% cl_location_selectivity; 
% 
% % movement stuffs
% cl_motionIndex; 
% 
% % rastermap/umap


% DMS vs PS 
