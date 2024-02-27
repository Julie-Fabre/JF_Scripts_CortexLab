
%% ~~ Naive ~~

%% Load data
load_type = 'naive';
loadVids = 0;
for iExperimentType = 1:6 %5:6% 1:6
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/naive_data', num2str(iExperimentType), '.mat'], '-struct', 'expData', '-v7.3');
end

close all;
load_type = 'taskGo';
loadVids = 0;
for iExperimentType = 2 %1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/gogogo_data', num2str(iExperimentType), '.mat'], '-struct', 'expData', '-v7.3');
end%80:error

close all;
load_type = 'taskNoGo';
loadVids = 0;
for iExperimentType = 2 %1:2
    [expData, session_data, regions] = cl_loadPerStimulusData(load_type, iExperimentType, loadVids);
    save(['/home/julie/Dropbox/MATLAB/goNogo_data', num2str(iExperimentType), '.mat'], '-struct', 'expData', '-v7.3');
end

%% Fig 1: visual responses in naïve mouse are sparsely present in striatum, more rare in GPe/SNr

%% -Example cells visual
% striatum 
close all;

unit=223

cl_plotExCell_psth('JF110', 2, 8, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'locations')


cl_plotExCell_psth('JF110', 1, 8, unit, 'stimOn_noMove', 2, ...
    [-0.2, 0.6], 0.001, 1, 1, 'spatialFreq')

cl_plotExCell_psth('JF110', 1, 8, unit, 'stimOn_noMove', 3, ...
    [-0.2, 0.6], 0.001, 1, 1, 'ori')


cl_plotExCell_psth('JF110', 3, 8, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'natImg')

% GPe
for unit = [61,57,83,47,42,40,8]
    cl_plotExCell_psth('JF109', 1,6, unit, 'stimOn_noMove', 2, ...
    [-0.2, 0.6], 0.001, 1, 1, 'spatialFreq')
end
%unit=8

cl_plotExCell_psth('JF109', 2, 6, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'locations')


cl_plotExCell_psth('JF109', 1,6, unit, 'stimOn_noMove', 2, ...
    [-0.2, 0.6], 0.001, 1, 1, 'spatialFreq')

cl_plotExCell_psth('JF109', 1, 6, unit, 'stimOn_noMove', 3, ...
    [-0.2, 0.6], 0.001, 1, 1, 'ori')


cl_plotExCell_psth('JF109', 3, 6, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'natImg')

% SNr 
% cl_plotExCell_psth(animal, experiment, iProbe, unit, align_type, group_type, ...
%    raster_window, psth_bin_size, plot_raster, plot_psth, color_type)
unit=30;
cl_plotExCell_psth('JF101', 1, 9, unit, 'stimOn_noMove', 2, ...
    [-0.2, 0.6], 0.001, 1, 1, 'spatialFreq')

cl_plotExCell_psth('JF101', 1, 9, unit, 'stimOn_noMove', 3, ...
    [-0.2, 0.6], 0.001, 1, 1, 'ori')

cl_plotExCell_psth('JF101', 2, 9, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'locations')

cl_plotExCell_psth('JF101', 3, 9, unit, 'stimOn_noMove', 1, ...
    [-0.2, 0.6], 0.001, 1, 1, 'natImg')

%% -Population cells visual
% redo 
cl_population_PSTH;

%% -Percentage cells visual
% using peak/trough 
cl_percentage_cells;


%% -location plots
%redo with peak 
pcells = false;
load_type = 'naive';
cl_location_PSTH_slice;

pcells = true;
cl_location_PSTH_slice;

%% allen atlas connectivity

%%  Fig 2: naïve striatal visual responses are selective to stimulus features
%% - Example cells sp freq.
%% - Example cells ori.
%% - Example cells location
%% - Example cells nat images

%% - selectivity index (simple and full)
regions = {'CP', 'GPe', 'SNr'};
cl_selectivity;

% Cell types
% fast responses in GPe/SNr ?

%% Fig 3: training in 3go task increases stimulus responses across the basal ganglia
%% - task performance 
%% - example cells
%% - percentage cells 
cl_fraction_cells_summary;

%% - location plot
pcells=false;
load_type='taskGo';
cl_location_PSTH_slice;

pcells=true;
cl_location_PSTH_slice;

pcells = false;
load_type = 'taskNaive';
cl_location_PSTH_slice;

pcells = true;
cl_location_PSTH_slice;


%% Fig 4: visual responses after 3go training are distinguishable in striatum, not GPe/SNr 
%% - example cells
%% - PSTHs
cl_plot_PSTHs(1, 0); %contra
%% - dot plot?
cl_stimId_gogogo;

%% - sel index? 

%% Fig 5: visual responses in go/no-go task are restricted to just “go” stimuli 
%% - percentage cells
%% - example cells 
%% - PSTHs
cl_plot_PSTHs(1, 0); %contra
%% - dot plots 
%% - drpime? sel? 




%% Supp.
% PSTHs

cl_plot_PSTHs(0, 1); %center


pcells=false;
load_type='taskNoGo';
cl_location_PSTH_slice;
pcells=true;
cl_location_PSTH_slice;
%% ~~ Suppl.: motion, movement ect. ; cell types ; delineating subregions ~~

%% Suppl: histo

%% Suppl: bombcell

%% others: RPE, ect.