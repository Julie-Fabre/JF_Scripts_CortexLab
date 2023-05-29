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
keep regions passive_data_per_cond

    responsive_cells = passive_data_per_cond.pvalue{1} < 0.05/4 |...
        passive_data_per_cond.pvalue{2} < 0.05/4 |...
        passive_data_per_cond.pvalue{3} < 0.05/4|...
        passive_data_per_cond.pvalue{4} < 0.05/4;

for iRegion=1:size(regions,2)
    responsive_cell_region(iRegion) = sum(responsive_cells' & passive_data_per_cond.unit_area ==iRegion) / ...
        sum(passive_data_per_cond.unit_area ==iRegion);
end

% increase vs decrease types 


figure();
plot()


