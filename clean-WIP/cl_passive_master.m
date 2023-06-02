%% passive_master script 

% do cells encode visual stimuli? 
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
passive_data_per_cond = cl_loadPerStimulusData('passive');
cl_selectivity;
%decoding? 
cl_location_selectivity; 

% movement stuffs
cl_motionIndex; 

% rastermap/umap


% DMS vs PS 
