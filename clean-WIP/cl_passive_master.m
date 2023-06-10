%% passive_master script 

% do cells encode visual stimuli? 
regions = {'CP', 'GPe', 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};

reload = 0; 
if reload 
    passive_data = cl_loadAverageData('passive'); % check we don't have any duplicates 
else
    passive_data = load();%passive.mat
end
cl_population_PSTH;
cl_location_PSTH;
cl_percentage_cells;

% cell classification 
cl_population_PSTH_celltype;
cl_celltype_example;
cl_GPe_celltype_playground; 

% cell selctivity/encoding of specific visual stimuli 
passive_data_gr = cl_loadPerStimulusData('passive', 1);
passive_data_lo = cl_loadPerStimulusData('passive', 2);
passive_data_nat1 = cl_loadPerStimulusData('passive', 3);
passive_data_nat2 = cl_loadPerStimulusData('passive', 6);
%choiceworld only 
cl_selectivity;
%decoding? 
cl_location_selectivity; 

% movement stuffs
cl_motionIndex; 

% rastermap/umap

% passive (for vs task) 
passive_data_cw1 = cl_loadPerStimulusData('passive', 4);
passive_data_cw2 = cl_loadPerStimulusData('passive', 5);

cl_plot_average_passive_cw; 
% DMS vs PS 
