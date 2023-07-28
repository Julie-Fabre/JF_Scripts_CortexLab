%% go no go guys 

bhvData = noGoWorld_behavior_debug({'JF093', 'JF091', 'JF104', 'JF105', 'JF107'});

for iMouse = 1:5
    bhvData(iMouse).stim_to_moveMean(bhvData(iMouse).stim_to_moveMean == 0) = NaN;
    rt(iMouse, :) = nanmean(bhvData(iMouse).stim_to_moveMean);
    
end

figure();
hold on;
mouseColors = {rgb('DarkOrchid'), rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black')}
for iMouse = 5:-1:1
plot(rt(iMouse, :), 'Color', mouseColors{iMouse});
end
xticks([1,2,3])
xticklabels({'stim 1', 'stim 2', 'stim 3'})
xlabel('stim type')
ylabel('mean reation time (s)')
makepretty
legend()

for iMouse = 1:5
    bhvData(iMouse).goLeft(bhvData(iMouse).goLeft == 0) = NaN;
    goLeft(iMouse, :) = nanmean(bhvData(iMouse).goLeft./bhvData(iMouse).nTrials);
    
end

figure();
hold on;
mouseColors = {rgb('DarkOrchid'), rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black')}
for iMouse = 5:-1:1
plot(goLeft(iMouse, :), 'Color', mouseColors{iMouse});
end
xticks([1,2,3])
xticklabels({'stim 1', 'stim 2', 'stim 3'})
xlabel('stim type')
ylabel('frac go')
makepretty
legend()



%% go go go guys 


bhvData = noGoWorld_behavior_debug({'JF097', 'JF096', 'JF106', 'JF099'});

for iMouse = 1:4
    %bhvData(iMouse).stim_to_moveMean(bhvData(iMouse).stim_to_moveMean == 0) = NaN;
    rt_gogogo(iMouse, :) = nanmean(bhvData(iMouse).stim_to_moveMean(end-11:end,:));
    
end

figure();
hold on;
mouseColors = {rgb('red'), rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black')}
for iMouse = 4:-1:1
plot(rt_gogogo(iMouse, :), 'Color', mouseColors{iMouse});
end
xticks([1,2,3])
xticklabels({'stim 1', 'stim 2', 'stim 3'})
xlabel('stim type')
ylabel('mean reaction time (s)')
makepretty
legend()
ylim([0.2, 0.9])

for iMouse = 1:4
    bhvData(iMouse).goLeft(bhvData(iMouse).goLeft == 0) = NaN;
    goLeft_gogogo(iMouse, :) = nanmean(bhvData(iMouse).goLeft./bhvData(iMouse).nTrials);
    
end

figure();
hold on;
mouseColors = {rgb('red'), rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black')}
for iMouse = 4:-1:1
plot(goLeft_gogogo(iMouse, :), 'Color', mouseColors{iMouse});
end
xticks([1,2,3])
xticklabels({'stim 1', 'stim 2', 'stim 3'})
xlabel('stim type')
ylabel('frac go')
makepretty
ylim([0 1])
legend()