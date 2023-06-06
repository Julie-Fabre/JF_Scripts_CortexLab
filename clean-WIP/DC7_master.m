%% DC7 _ master 


%% passive stuff
% averages 
passive_data = cl_loadAverageData('passive');
cl_population_PSTH;
cl_location_PSTH;
cl_percentage_cells;
% example cells: 


% passive (for vs task) 
passive_data_cw1 = cl_loadPerStimulusData('passive', 4);
passive_data_cw2 = cl_loadPerStimulusData('passive', 5);

cl_plot_average_passive_cw; 
%% task stuff 
task_data = cl_loadPerStimulusData('task', 2); 
cl_plot_average_task; 