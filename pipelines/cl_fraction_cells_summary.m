close all;
keepVis = 0;
passive = 1;
goNogo = 0;
cl_fraction_cells;

save('increaseFR_session_fraction_passive.mat', 'increaseFR_session_fraction')

passive = 0;
goNogo = 0;
cl_fraction_cells;

save('increaseFR_session_fraction_gogogo.mat', 'increaseFR_session_fraction')


passive = 0;
goNogo = 1;
cl_fraction_cells;

save('increaseFR_session_fraction_gonogo.mat', 'increaseFR_session_fraction')

%%
fraction_passive = load('increaseFR_session_fraction_passive.mat');
fraction_gonogo = load('increaseFR_session_fraction_gonogo.mat');
fraction_gogogo = load('increaseFR_session_fraction_gogogo.mat');
%(iRegion,iSession, iPair)
cols = lines(3); 
figure();
regions={'Str', 'GPe', 'SNr'};
for iRegion=1:3
subplot(1,3,iRegion); hold on;

swarmchart(ones(size(fraction_passive.increaseFR_session_fraction,2),1), ...
    fraction_passive.increaseFR_session_fraction(iRegion, :), 20, cols(1,:));
swarmchart(ones(size(fraction_gonogo.increaseFR_session_fraction,2),1).*3, ...
    fraction_gonogo.increaseFR_session_fraction(iRegion, :), 20, cols(3,:));
swarmchart(ones(size(fraction_gogogo.increaseFR_session_fraction,2),1).*2, ...
    fraction_gogogo.increaseFR_session_fraction(iRegion, :), 20, cols(2,:));
ylim([0, 0.35])
title(regions{iRegion})

end
subplot(1,3,1);
ylabel('% visual cells')
prettify_plot;

figure();
for iRegion=1:3
subplot(1,3,iRegion); hold on;

    bar([1,2,3],[nanmean(fraction_passive.increaseFR_session_fraction(iRegion, :, 1)), ...
        nanmean(fraction_gogogo.increaseFR_session_fraction(iRegion, :, 1)), nanmean(fraction_gonogo.increaseFR_session_fraction(iRegion, :, 1))])
%ylim([0, 0.45])
end
subplot(1,3,1);
ylabel(['mean % cells' newline 'd-prime > 0.5 per session'])
xlabel('Pair # ')
prettify_plot

for iRegion =1:3

ranksum(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gogogo.increaseFR_session_fraction(iRegion, :))
ranksum(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :))
ranksum(fraction_gogogo.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :))

%ttest2(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gogogo.increaseFR_session_fraction(iRegion, :))
%ttest2(fraction_passive.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :))
%ttest2(fraction_gogogo.increaseFR_session_fraction(iRegion, :), fraction_gonogo.increaseFR_session_fraction(iRegion, :))



end
