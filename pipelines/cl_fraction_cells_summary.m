clear all;
close all;
keepVis = 0;
passive = 1;
goNogo = 0;
cl_fraction_cells;

save('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_passive.mat', 'increaseFR_session_fraction')

keep keepVis
passive = 0;
goNogo = 0;
cl_fraction_cells;

save('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_gogogo.mat', 'increaseFR_session_fraction')

keep keepVis
passive = 0;
goNogo = 1;
cl_fraction_cells;

save('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_gonogo.mat', 'increaseFR_session_fraction')

%%
fraction_passive = load('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_passive.mat');
fraction_gonogo = load('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_gonogo.mat');
fraction_gogogo = load('/home/julie/Dropbox/MATLAB/increaseFR_session_fraction_gogogo.mat');
%(iRegion,iSession, iPair)
cols = brewermap(3, 'Set1');
figure();
regions = {'Str', 'GPe', 'SNr'};
for iRegion = 1:3
    subplot(1, 3, iRegion);
    hold on;

    violinplot([fraction_passive.increaseFR_session_fraction(iRegion, :) .* 100, ...
        fraction_gonogo.increaseFR_session_fraction(iRegion, :) .* 100, ...
        fraction_gogogo.increaseFR_session_fraction(iRegion, :) .* 100], ...
        [ones(size(fraction_passive.increaseFR_session_fraction, 2), 1); ...
        ones(size(fraction_gonogo.increaseFR_session_fraction, 2), 1) .* 3; ...
        ones(size(fraction_gogogo.increaseFR_session_fraction, 2), 1) .* 2], ...
        'ViolinColor', cols, 'ShowMedian', true, 'MedianColor', [1,1,1],...
        'MedianMarkerSize', 45, 'ShowBox', false, 'ShowWhiskers', true);

    %ylim([-0.1, 0.35])
    title(regions{iRegion})
    xticks([1, 2, 3])
    xticklabels({'Naive', 'Gogogo', 'GoNogo'})
    ylabel('% visual cells');
end
prettify_plot('YLimits', 'all');

for iRegion = 1:3
    subplot(1, 3, iRegion);
    hold on;
    [~, pval_12] = ttest(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gogogo.increaseFR_session_fraction(iRegion, :));
    [~, pval_13] = ttest(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :));
    [~, pval_23] = ttest(fraction_gogogo.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :));
    disp([pval_12, pval_13, pval_23])

    prettify_pvalues(gca, [1, 1, 2], [2, 3, 3], [pval_12, pval_13, pval_23],'PlotNonSignif',false)
end

