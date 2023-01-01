%% JF _ overview of meta functions 

% loading data
allDataStruct = JF_loadAllData(animals, region, protocol, useQualityMetrics);

% kilosorting, preprocessing, decompressing 

% quality metrics

% histology

% PSTHs
[curr_smoothed_psth, curr_psth, raster_x, raster_y, rasterColor, alignedVector] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group, sort_by, color_by, plot_me, causal_smoothing, figN);
allData_singleCellPSTH = JF_singleCellPSTH_allData(allDataStruct,trialGroups, alignTo, raster_window, psth_bin_size);
JF_singleCellPSTH_flip_allData(allDataStruct,trialGroups, alignTo, raster_window, psth_bin_size)

% behavior 
bhvOut = JF_noGoWorld_behavior(animalsAll);
% false alarm, hit, correct rates

% pupil/movement tracking

