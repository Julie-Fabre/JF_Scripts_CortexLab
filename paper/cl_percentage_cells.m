
%% cl_percentage_cells
%
% from : https://www.nature.com/articles/s41467-021-25436-3#Sec10
% To characterize a neuronfs responsiveness to a given visual stimulus, we
% calculated its “reliability” on every session56, defined by the Pearson
% correlation coefficient (CC) of the session-averaged activity of two
% random halves of trials, iterated many times, and averaged (see “Methods”).
% Reliability for both stimuli followed a skewed distribution, with a higher
% number of neurons responding reliably to the MOV stimulus (Fig. 1d). To
% determine which neurons were visually responsive we tested the actual
% reliability of the response to each stimulus against a null distribution
% of reliability calculated with circularly shifted data
% A neuron%s responsiveness to a stimulus was determined based on a collective
% measure of the reliability of the neurons in a given field using time-shuffled
% data. First, a neuron s activity on each trial was circularly shuffled by a
% random amount. Next, a reliability value was calculated using this shuffled
% data. This was repeated 1000 times to yield a distribution of reliability
% values, and the 99th percentile of this distribution was stored. This 99th
% percentile threshold was found for every neuron. If a neuron’s actual average
% reliability across sessions was statistically greater than the average of
% these 99th percentile values (two-tailed one-sample t-test), it was classified
% as responsive to the stimulus.

% just do a sign rank test (in cl_loadPerStimulusData)(positive for any of
% conditions w/ holm-bonferroni test)
passive_data = load('/home/julie/Dropbox/MATLAB/naive_data1.mat');
keep regions task_data passive_data
responsive_cells = passive_data.pvalue < 0.05 & ...
    passive_data.pvalue_shuffled == 1;
zscore_psth = (passive_data.psth - nanmean(passive_data.psth(:, 1:50), 2)) ./ ...
    nanstd(passive_data.psth(:, 1:50), [], 2);
keep_these = sum(isnan(zscore_psth), 2) < 100;
for iRegion = 1:size(regions, 2)
    responsive_cell_region(iRegion) = sum(responsive_cells' & passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & keep_these) / ...
        sum(passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & ...
        keep_these);

    responsive_cell_region_up(iRegion) = sum(responsive_cells' & passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & keep_these & nanmean(zscore_psth(:, 55:65), 2) > 0) / ...
        sum(passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & ...
        keep_these);

    responsive_cell_region_down(iRegion) = sum(responsive_cells' & passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & keep_these & nanmean(zscore_psth(:, 55:65), 2) < 0) / ...
        sum(passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & ...
        keep_these);


end

cl_plottingSettings;
colorMtx = bc_colors(4);
figure();
b = bar(1:3, [responsive_cell_region([1, 2, 5])], 0.2, 'facecolor', 'flat');
b.CData = [regionColors{1}; regionColors{2}; regionColors{3}];
xticks([1:3])
xticklabels({'STR', 'GPe', 'SNr'})
makepretty;
ylabel(['fraction of', newline, 'visual cells'])
xlim([0.5, 3.5])

iRegion = 1;
cl_plottingSettings;
colorMtx = bc_colors(4);
figure();
b = barh(1:2, [responsive_cell_region_up(iRegion), responsive_cell_region_down(iRegion)], 0.2, 'facecolor', 'flat');
b.CData = [regionColors{iRegion}; regionColors{iRegion}];
yticks([1:2])
yticklabels({'decrease', 'increase'})
makepretty;
ylabel(['fraction of', newline, 'visual cells'])
ylim([0.5, 2.5])
xlim([0, 0.076])


responsive_cells_task = task_data.pvalue{2} < 0.05 & ...
    task_data.pvalue_shuffled_005{2} == 1;
zscore_psth = squeeze((task_data.psth{2}(:, 1, 1, :))-nanmean(squeeze(task_data.psth{2}(:, 1, 1, 1:50)), 2)) ./ ...
    nanstd(squeeze(task_data.psth{2}(:, 1, 1, 1:50)), [], 2);
keep_these = sum(isnan(zscore_psth), 2) < 100;
for iRegion = 1:size(regions, 2)
    responsive_cell_region_task(iRegion) = 2 * sum(responsive_cells_task' & task_data.unit_area == iRegion & ...
        (task_data.unitType' == 1 | task_data.unitType' == 2)) / ...
        sum(task_data.unit_area == iRegion & ...
        (task_data.unitType' == 1 | task_data.unitType' == 2) & ...
        keep_these);
end

% increase vs decrease types

cl_plottingSettings;
colorMtx = bc_colors(4);
figure();
b = barh(2:-1:1, [responsive_cell_region(1), responsive_cell_region_task(1)], 0.2, 'facecolor', 'flat');
b.CData = [colorMtx(4, :); colorMtx(3, :)];
yticks([1:2])
yticklabels({'pre-learning', 'post-learning'})
makepretty;
xlabel(['fraction of', newline, 'visual cells'])
