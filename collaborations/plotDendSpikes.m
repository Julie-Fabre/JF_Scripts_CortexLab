%% data from Giulia 
% 
load('/home/julie/Dropbox/Data/Shared/Giulia/map_unit117.mat')
load('/home/julie/Dropbox/Data/Shared/Giulia/map_unit134.mat')

waveform_time = 1e3 * ((0:size(map_unit117, 2) - 1) / 30000);


[~, somaLocation] = min(min(map_unit117,[],2));
%maxVal = max(max(abs(map_unit117)));
maxVal = 10;
waveformLimits = [find(waveform_time > 1, 1, 'first'), find(waveform_time > 13, 1, 'first')];
%siteSpacing = 
figure();
hold on;
for iChannel = somaLocation - 10:somaLocation + 50
    if iChannel == somaLocation
        plot(waveform_time(waveformLimits(1):waveformLimits(2)) - waveform_time(waveformLimits(1)), ...
            map_unit117(iChannel, waveformLimits(1):waveformLimits(2)) + maxVal * ...
            (iChannel - somaLocation), 'Color', rgb('Cyan'), 'LineWidth',4)

    else
        plot(waveform_time(waveformLimits(1):waveformLimits(2)) - waveform_time(waveformLimits(1)), ...
            map_unit117(iChannel, waveformLimits(1):waveformLimits(2)) + maxVal * ...
            (iChannel - somaLocation), 'Color', 'w')
    end
end
yticks([0, 200, 400])
yticklabels({'0', '400', '800'})
ylabel('distance from soma')
xlabel('time (ms)') 
prettify_plot('','','k');
ylim([-110, 338])
xlim([4, 8])
prettify_addScaleBars([], [], [], [], [], 'ms', 'um');


figure();
colormap(brewermap([], '*RdBu'));
imagesc(flipud(map_unit134))
clim([-20, 20])



%% plot noise / good 
bc_qualityMetricsPipeline_JF('JF093','2023-03-06',1,[],1,[],1,1,1)

noiseUnits = find(unitType == 0);
channelPositions_x = repmat([43,11,59,27], 1, 96);
channelPositions_y = sort(repmat(0:20:3820, 1, 2));
waveform_time = 1e3 * ((0:size(rawWaveformsFull, 3) - 1) / 30000);
scaleF = [repmat(1, 24, 1); 2];

figure('Color', 'k');
clf;
noiseUnits = [2,21,25,11];
for idx = 1:4
iNoiseUnit = noiseUnits(idx);
subplot(1,4,idx)
hold on;
set(gca, 'color', 'k');
for iChannel = 1:size(channelPositions,1)
plot(channelPositions_x(iChannel)./6 + waveform_time, ...
squeeze(rawWaveformsFull(iNoiseUnit,iChannel,:)* scaleF(iNoiseUnit)) + ...
(channelPositions_y(iChannel) ), 'Color', 'w')
end
xlim([1.8, 12.6])
ylim([-50, 3600])
end
prettify_plot('','','k');
for idx = 1:4
subplot(1,4,idx)
prettify_addScaleBars([], [], [], [], [], 'ms', 'um');
end