%% passive_master script 

% do cells encode visual stimuli? 
cl_loadAverageData;
cl_population_PSTH;
cl_location_PSTH;
cl_percentage_cells;

% cell classification 
cl_population_PSTH_celltype;

% cell selctivity/encoding of specific visual stimuli 
cl_loadPerStimulusData;
cl_selectivity;
%decoding? 
cl_location_selectivity; 

% movement stuffs
cl_motionIndex; 

% rastermap/umap


% DMS vs PS 