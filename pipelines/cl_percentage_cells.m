% from https://www.nature.com/articles/s41467-021-25436-3
% A neuronfs responsiveness to a stimulus was determined based on a collective 
% measure of the reliability of the neurons in a given field using time-shuffled 
% data. First, a neuron’s activity on each trial was circularly shuffled by 
% a random amount. Next, a reliability value was calculated using this shuffled 
% data. This was repeated 1000 times to yield a distribution of reliability 
% values, and the 99th percentile of this distribution was stored. This 99th 
% percentile threshold was found for every neuron. If a neuron’s actual 
% average reliability across sessions was statistically greater than the 
% average of these 99th percentile values (two-tailed one-sample t-test), 
% it was classified as responsive to the stimulus. 

passive = 0;
goNogo = 0;
keepVis = 1;
cl_fraction_cells;


cols = brewermap(3, 'Set1');
figure();
regions = {'Str', 'GPe', 'SNr'};
%for iRegion = 1:3
  %  subplot(1, 3, iRegion);
    hold on;

    violinplot([increaseFR_session_fraction(1, :) .* 100, ...
        increaseFR_session_fraction(2, :) .* 100, ...
        increaseFR_session_fraction(3, :) .* 100], ...
        [ones(size(increaseFR_session_fraction, 2), 1); ...
        ones(size(increaseFR_session_fraction, 2), 1) .* 3; ...
        ones(size(increaseFR_session_fraction, 2), 1) .* 2], ...
        'ViolinColor', cols, 'ShowMedian', true, 'MedianColor', [1,1,1],...
        'MedianMarkerSize', 45, 'ShowBox', false, 'ShowWhiskers', true);

    %ylim([-0.1, 0.35])
    title(regions{iRegion})
    xticks([1, 2, 3])
    xticklabels({'CP', 'GPe', 'SNr'})
    ylabel('% visual cells');
%end
prettify_plot;
