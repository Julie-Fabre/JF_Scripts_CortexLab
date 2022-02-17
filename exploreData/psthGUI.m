

function psthGUI(task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
    task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
    passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
    passive_wheel_types, passive_trial_conditions)

%% set up dynamic figure
psthGuiHandle = figure('color', 'w');
set(psthGuiHandle, 'KeyPressFcn', @KeyPressCb);

%% initial conditions
iCluster = 1;
uniqueTemps = unique(ephysData.spike_templates);

%% plot initial conditions
initializePlot(psthGuiHandle)
updateUnit(psthGuiHandle,task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
    task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
    passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
    passive_wheel_types, passive_trial_conditions);

%% change on keypress
    function KeyPressCb(~, evnt)
        %fprintf('key pressed: %s\n', evnt.Key);
        if strcmpi(evnt.Key, 'rightarrow')
            iCluster = iCluster + 1;
            updateUnit(psthGuiHandle,task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
                task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
                passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
                passive_wheel_types, passive_trial_conditions);
        elseif strcmpi(evnt.Key, 'leftarrow')
            iCluster = iCluster - 1;
            updateUnit(psthGuiHandle,task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
                task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
                passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
                passive_wheel_types, passive_trial_conditions);
        end
    end
end

function initializePlot(psthGuiHandle)

%% main title

mainTitle = sgtitle('');

%% initialize task stim PSTH 

subplot(6, 2, [1, 3]);
hold on;
rasterTaskStimDots = scatter(NaN,NaN,4,'k','filled');
rasterTaskStimLine = line([NaN, NaN], [NaN, NaN])
rasterTaskMoveLine
rasterTaskRewaLine 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')

%% initialize task stim PSTH 
subplot(6, 2, 5);
hold on;
psthTaskLines = scatter(NaN,NaN,4,'k','filled');
psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')

%% initialize task spontaneous move PSTH
subplot(6, 2, [2, 4])
hold on;
psthTaskLines = scatter(NaN,NaN,4,'k','filled');
psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')

%% initialize 
subplot(2, 2, 3)
hold on;
rawWaveformLines = arrayfun(@(x) plot(NaN, NaN, 'linewidth', 2, 'color', 'k'), 1:max_n_channels_plot);
maxRawWaveformLines = arrayfun(@(x) plot(nan(82, 1), nan(82, 1), 'linewidth', 2, 'color', 'b'), 1);
set(gca, 'YDir', 'reverse')
xlabel('Position+Time');
ylabel('Position');
rawTitle = title('');
rawLegend = legend([maxRawWaveformLines], {''});

%% initialize ACG
subplot(2, 2, 4)

hold on;
acgBar = arrayfun(@(x) bar(0:0.1:25, nan(251, 1)), 1);
acgRefLine = line([NaN, NaN], [NaN, NaN], 'Color', 'r', 'linewidth', 1.2);
acgAsyLine = line([NaN, NaN], [NaN, NaN], 'Color', 'r', 'linewidth', 1.2);
xlabel('time (ms)');
ylabel('sp/s');
acgTitle = title('');


%% save all handles
guiData = struct;
% main title
guiData.mainTitle = mainTitle;
% location plot
guiData.unitDots = unitDots;

% upload guiData
guidata(unitQualityGuiHandle, guiData);
end

function updateUnit(psthGuiHandle,task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
    task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
    passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
    passive_wheel_types, passive_trial_conditions)

%% Get guidata
guiData = guidata(unitQualityGuiHandle);
thisUnit = uniqueTemps(iCluster);

%% main title
if unitType(iCluster) == 1
    set(guiData.mainTitle, 'String', ['Unit ', num2str(iCluster), ', single unit'], 'Color', [0, .5, 0]);
elseif unitType(iCluster) == 0
    set(guiData.mainTitle, 'String', ['Unit ', num2str(iCluster), ', noise/non-somatic'], 'Color', [1, 0, 0]);
elseif unitType(iCluster) == 2
    set(guiData.mainTitle, 'String', ['Unit ', num2str(iCluster), ', multi-unit'], 'Color', [0.29, 0, 0.51]);
end

%% plot 1: update curr unit location
set(guiData.currUnitDots, 'XData', guiData.norm_spike_n(thisUnit), 'YData', ephysData.channel_positions(qMetrics.maxChannels(thisUnit), 2), 'CData', guiData.unitCmap(iCluster, :))

%% plot 2: update unit template waveform and detected peaks
% guiData.templateWaveformLines = templateWaveformLines;
%     guiData.maxTemplateWaveformLines = maxTemplateWaveformLines;
%     guiData.tempTitle = tempTitle;
%     guiData.tempLegend = tempLegend;

maxChan = qMetrics.maxChannels(thisUnit);
maxXC = ephysData.channel_positions(maxChan, 1);
maxYC = ephysData.channel_positions(maxChan, 2);
chanDistances = ((ephysData.channel_positions(:, 1) - maxXC).^2 ...
    +(ephysData.channel_positions(:, 2) - maxYC).^2).^0.5;
chansToPlot = find(chanDistances < 100);
vals =[];
for iChanToPlot = 1:min(20, size(chansToPlot, 1))
    vals(iChanToPlot) = max(abs(squeeze(ephysData.templates(thisUnit, :, chansToPlot(iChanToPlot)))));
    if maxChan == chansToPlot(iChanToPlot)
        set(guiData.maxTemplateWaveformLines, 'XData', (ephysData.waveform_t + (ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(ephysData.templates(thisUnit, :, chansToPlot(iChanToPlot)))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) ./ 100));
        hold on;
        set(guiData.peaks, 'XData', (ephysData.waveform_t(qMetrics.peakLocs{iCluster}) ...
            +(ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(ephysData.templates(thisUnit, qMetrics.peakLocs{iCluster}, chansToPlot(iChanToPlot)))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) ./ 100));

        set(guiData.troughs, 'XData', (ephysData.waveform_t(qMetrics.troughLocs{iCluster}) ...
            +(ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(ephysData.templates(thisUnit, qMetrics.troughLocs{iCluster}, chansToPlot(iChanToPlot)))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) ./ 100));
        set(guiData.templateWaveformLines(iChanToPlot), 'XData', nan(82, 1), ...
            'YData', nan(82, 1));

    else
        set(guiData.templateWaveformLines(iChanToPlot), 'XData', (ephysData.waveform_t + (ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(ephysData.templates(thisUnit, :, chansToPlot(iChanToPlot)))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) ./ 100));
    end
end
%  [nPeaks, nTroughs, somatic, peakLocs, troughLocs] = bc_troughsPeaks(ephysData.templates(thisUnit, :, qMetrics.maxChannels(thisUnit)), ...
%          param.ephys_sample_rate, 1);

% ff=find(chansToPlot == maxChan);
% disp(vals(ff))
% if ff>4 && ff<min(20, size(chansToPlot, 1))-4
%     disp(nanmean(vals([ff-2, ff-1, ff+1, ff+2])))
%     disp(nanmean(vals([ff-4, ff-3, ff+4, ff+3])))
% elseif ff<min(20, size(chansToPlot, 1))-4
%     disp(nanmean(vals([ff+1, ff+2])))
%     disp(nanmean(vals([ff+4, ff+3])))
% elseif ff>4
%     disp(nanmean(vals([ff-2, ff-1])))
%     disp(nanmean(vals([ff-4, ff-3])))
% end

% 
% figure();    
% scatter3(vals, ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot, 1))),1), ...
%     ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot,1))),2))

% X_ave=mean([vals', ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot, 1))),1), ...
%     ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot,1))),2)],1);            
% % mean; line of best fit will pass through this point  
% dX=bsxfun(@minus,[vals', ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot, 1))),1), ...
%     ephysData.channel_positions(chansToPlot(1:min(12, size(chansToPlot,1))),2)],X_ave);  % residuals
% C=(dX'*dX)/(numel(vals)-1);           % variance-covariance matrix of X
% [R,D]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
% D=diag(D);
% R2=D(1)/sum(D);

if qMetrics.nPeaks(iCluster) > param.maxNPeaks || qMetrics.nTroughs(iCluster) > param.maxNTroughs
    if qMetrics.somatic(iCluster) == 0
        set(guiData.tempTitle, 'String', ['\fontsize{9}Template waveform: {\color[rgb]{1 0 0}# detected peaks/troughs, ', ...
            '\color[rgb]{1 0 0}is somatic \color{red}}'])
    else
        set(guiData.tempTitle, 'String', ['\fontsize{9}Template waveform: {\color[rgb]{1 0 0}# detected peaks/troughs, ', ...
            '\color[rgb]{0 .5 0}is somatic \color{red}}'])
    end
else
    if qMetrics.somatic(iCluster) == 0
        set(guiData.tempTitle, 'String', ['\fontsize{9}Template waveform: {\color[rgb]{0 .5 0}# detected peaks/troughs, ', ...
            '\color[rgb]{1 0 0}is somatic \color{red}}'])
    else
        set(guiData.tempTitle, 'String', ['\fontsize{9}Template waveform: {\color[rgb]{0 .5 0}# detected peaks/troughs, ', ...
            '\color[rgb]{0 .5 0}is somatic \color{red}}'])
    end

end
set(guiData.tempLegend, 'String', {['is somatic =', num2str(qMetrics.somatic(iCluster))], ...
    [num2str(qMetrics.nPeaks(iCluster)), ' peak(s)'], [num2str(qMetrics.nTroughs(iCluster)), ' trough(s)']})

%% plot 3: plot unit mean raw waveform (and individual traces)

for iChanToPlot = 1:min(20, size(chansToPlot, 1))
    if maxChan == chansToPlot(iChanToPlot)
        set(guiData.maxRawWaveformLines, 'XData', (ephysData.waveform_t + (ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(qMetrics.rawWaveforms(iCluster).spkMapMean(chansToPlot(iChanToPlot), :))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) * 10));
        set(guiData.rawWaveformLines(iChanToPlot), 'XData', nan(82, 1), ...
            'YData', nan(82, 1));

    else
        set(guiData.rawWaveformLines(iChanToPlot), 'XData', (ephysData.waveform_t + (ephysData.channel_positions(chansToPlot(iChanToPlot), 1) - 11) / 10), ...
            'YData', -squeeze(qMetrics.rawWaveforms(iCluster).spkMapMean(chansToPlot(iChanToPlot), :))'+ ...
            (ephysData.channel_positions(chansToPlot(iChanToPlot), 2) * 10));
    end
end
set(guiData.rawLegend, 'String', ['Amplitude =', num2str(qMetrics.rawAmplitude(iCluster)), 'uV'])
if qMetrics.rawAmplitude(iCluster) < param.minAmplitude
    set(guiData.rawTitle, 'String', '\color[rgb]{1 0 1}Mean raw waveform');
else
    set(guiData.rawTitle, 'String', '\color[rgb]{0 .5 0}Mean raw waveform');
end

%% 4. plot unit ACG

theseSpikeTimes = ephysData.spike_times_timeline(ephysData.spike_templates == thisUnit);

[ccg, ccg_t] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
    ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', 0.001, 'duration', 0.5, 'norm', 'rate'); %function

set(guiData.acgBar, 'XData', ccg_t(250:501)*1000, 'YData', squeeze(ccg(250:501, 1, 1)));
set(guiData.acgRefLine, 'XData', [2, 2], 'YData', [0, max(ccg(:, 1, 1))])
[ccg2, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
    ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', 0.1, 'duration', 10, 'norm', 'rate'); %function
asymptoteLine = nanmean(ccg2(end-100:end));
set(guiData.acgAsyLine, 'XData', [0, 250], 'YData', [asymptoteLine, asymptoteLine])

if qMetrics.Fp(iCluster) > param.maxRPVviolations
    set(guiData.acgTitle, 'String', '\color[rgb]{1 0 1}ACG');
else
    set(guiData.acgTitle, 'String', '\color[rgb]{0 .5 0}ACG');
end

%% 5. plot unit ISI (with refractory period and asymptote lines)

theseISI = diff(theseSpikeTimes);
theseISIclean = theseISI(theseISI >= param.tauC); % removed duplicate spikes
theseOffendingSpikes = find(theseISIclean < (2/1000)); 
theseOffendingSpikes = [theseOffendingSpikes; theseOffendingSpikes-1];
[isiProba, edgesISI] = histcounts(theseISIclean*1000, [0:0.5:50]);

set(guiData.isiBar, 'XData', edgesISI(1:end-1)+mean(diff(edgesISI)), 'YData', isiProba); %Check FR
set(guiData.isiRefLine, 'XData', [2, 2], 'YData', [0, max(isiProba)])

if qMetrics.Fp(iCluster) > param.maxRPVviolations
    set(guiData.isiTitle, 'String', '\color[rgb]{1 0 1}ISI');
else
    set(guiData.isiTitle, 'String', '\color[rgb]{0 .5 0}ISI');
end
set(guiData.isiLegend, 'String', [num2str(qMetrics.Fp(iCluster)), ' % r.p.v.'])

%% 6. plot isolation distance
if param.computeDistanceMetrics
set(guiData.currIsoD, 'XData', qMetrics.Xplot{iCluster}(:, 1), 'YData', qMetrics.Xplot{iCluster}(:, 2))
set(guiData.rpvIsoD, 'XData', qMetrics.Xplot{iCluster}(theseOffendingSpikes, 1), 'YData', qMetrics.Xplot{iCluster}(theseOffendingSpikes, 2))
set(guiData.otherIsoD, 'XData', qMetrics.Yplot{iCluster}(:, 1), 'YData', qMetrics.Yplot{iCluster}(:, 2), 'CData', qMetrics.d2_mahal{iCluster})
end
%% 7. (optional) plot raster

%% 10. plot ampli fit

    
set(guiData.ampliBins, 'XData', qMetrics.ampliBinCenters{iCluster}, 'YData', qMetrics.ampliBinCounts{iCluster});

set(guiData.ampliFit, 'XData', qMetrics.ampliFit{iCluster}, 'YData', qMetrics.ampliBinCenters{iCluster})
if qMetrics.percSpikesMissing(iCluster) > param.maxPercSpikesMissing
    set(guiData.ampliFitTitle, 'String', '\color[rgb]{1 0 1}% spikes missing');
else
    set(guiData.ampliFitTitle, 'String', '\color[rgb]{0 .5 0}% spikes missing');
end
set(guiData.ampliFitLegend, 'String', {[num2str(qMetrics.percSpikesMissing(iCluster)), ' % spikes missing'], 'rpv spikes'})
set(guiData.ampliFitAx, 'YLim', [min(qMetrics.ampliBinCenters{iCluster}), max(qMetrics.ampliBinCenters{iCluster})])

%% 9. plot template amplitudes and mean f.r. over recording (QQ: add experiment time epochs?)

ephysData.recordingDuration = (max(ephysData.spike_times_timeline) - min(ephysData.spike_times_timeline));
theseAmplis = ephysData.template_amplitudes(ephysData.spike_templates == thisUnit);

% for debugging if wierd amplitude fit results: bc_percSpikesMissing(theseAmplis, theseSpikeTimes, [min(theseSpikeTimes), max(theseSpikeTimes)], 1);

set(guiData.tempAmpli, 'XData', theseSpikeTimes, 'YData', theseAmplis)
set(guiData.rpvAmpli, 'XData', theseSpikeTimes(theseOffendingSpikes), 'YData', theseAmplis(theseOffendingSpikes))
currTimes = theseSpikeTimes(theseSpikeTimes >= theseSpikeTimes(iChunk)-0.1 & theseSpikeTimes <= theseSpikeTimes(iChunk)+0.1);
currAmplis = theseAmplis(theseSpikeTimes >= theseSpikeTimes(iChunk)-0.1 & theseSpikeTimes <= theseSpikeTimes(iChunk)+0.1);
set(guiData.currTempAmpli, 'XData', currTimes, 'YData', currAmplis);
set(guiData.ampliAx.YAxis(1), 'Limits', [0, round(max(theseAmplis))])

binSize = 20;
timeBins = 0:binSize:ceil(ephysData.spike_times(end)/ephysData.ephys_sample_rate);
[n, x] = hist(theseSpikeTimes, timeBins);
n = n ./ binSize;

set(guiData.spikeFR, 'XData', x, 'YData', n);
set(guiData.ampliAx.YAxis(2), 'Limits', [0, 2 * ceil(max(n))])


if qMetrics.nSpikes(iCluster) > param.minNumSpikes
    set(guiData.ampliTitle, 'String', '\color[rgb]{0 .5 0}Spikes');
else
    set(guiData.ampliTitle, 'String', '\color[rgb]{1 0 1}Spikes');
end
set(guiData.ampliLegend, 'String', {['# spikes = ', num2str(qMetrics.nSpikes(iCluster))], 'rpv spikes'})

%% 8. plot raw data
if plotRaw
    plotSubRaw(guiData.rawPlotH, guiData.rawPlotLines, guiData.rawSpikeLines, memMapData, ephysData, iCluster, uniqueTemps, iChunk);
end
end


