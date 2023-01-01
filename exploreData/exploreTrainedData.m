%% exploring trained data 
% no JF070 striatum ?? 
% STR: JF070 21 (bottom after 3200 - posterior striatum ++) ? ; bottom 31 after 3600 top  1 1
% (until 1500)- after maybe went through white matter trackt?  ? 
% JF071 (/!\ azimuth passive) 
% JF078
% JF067 (/!\ azimuth passive)
% GPe : JF070 day 2,3,4 site 1 (?) 
% SNr : JF070

%% striatum, chaque bonne cellule, baseline go/no go action 
allDataStruct = JF_loadAllData({'JF070'}, 'CP', 'stage', 1);
raster_window = [-0.5, 2];
psth_bin_size = 0.01;
allData_singleCellPSTH = JF_singleCellPSTH_allData(allDataStruct,'stim_id', 'stim', raster_window, psth_bin_size);

% plot baseline in go vs no go trials depending on stim type 
baseline_goStim_goAction = nanmean(allData_singleCellPSTH.psth(:,1,56:71),3)./1.5;
baseline_goStim_noGoAction = nanmean(allData_singleCellPSTH.psth(:,4,56:71),3)./1.5;
baseline_noGoStim_goAction = nanmean(allData_singleCellPSTH.psth(:,3,56:71),3)./1.5;
baseline_noGoStim_noGoAction = nanmean(allData_singleCellPSTH.psth(:,5,56:71),3)./1.5;
figure('Color', 'white');
scatterHistDiff(baseline_goStim_goAction(baseline_goStim_goAction<=30), baseline_goStim_noGoAction(baseline_goStim_noGoAction<=30), [], [], 'blue',1)
xlim([0 30])
ylim([0 30])
xlabel(['baseline firing rate in' newline ...
    'go stimulus, go action trials'])
ylabel(['baseline firing rate in' newline ...
    'go stimulus, no go action trials'])

figure('Color', 'white');
scatterHistDiff(baseline_noGoStim_goAction(baseline_noGoStim_goAction<=30), baseline_noGoStim_noGoAction(baseline_noGoStim_noGoAction<=30), [], [], 'blue',1)
xlim([0 30])
ylim([0 30])
xlabel(['baseline firing rate in' newline ...
    'no go stimulus, go action trials'])
ylabel(['baseline firing rate in' newline ...
    'no go stimulus, no go action trials'])

JF_singleCellPSTH_flip_allData(allDataStruct,'stim_id', 'stim', raster_window, psth_bin_size);

% -> FSIs. show example
% -> MSNs. show example 

% single trial population raster 

% correlation no go rate with pupil/whisking 

% example cell 
[curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group, sort_by, color_by, plot_me, causal_smoothing)

JF_singleTrialPSTH_allData

