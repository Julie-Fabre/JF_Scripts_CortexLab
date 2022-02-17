
%% load data 
animals = {'AP100','AP101','AP104','AP105','AP106'};
passiveProtocol = 'AP_lcrGratingPassive';
bhvProtocol = 'AP_stimWheelRight';

for iAnimal = 1:size(animals,2)
animal = animals{iAnimal};
experiments = AP_find_experimentsJF(animal, bhvProtocol, true);
experiments = experiments([experiments.ephys]);

for iDay = 1:size(experiments,1)
day = experiments(iDay).day;
experiment = experiments(iDay).experiment;

ephysPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys');
[spikeTimes, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, usedChannels] = bc_loadEphysData(ephysPath);
ephysap_path = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap');
ephysDirPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_dir');
savePath = fullfile(ephysDirPath, 'qMetrics'); 

%% run qmetrics 
param = struct;
param.plotThis = 0;
param.plotGlobal =1;
% refractory period parameters
param.tauR = 0.0020; %refractory period time (s)
param.tauC = 0.0001; %censored period time (s)
param.maxRPVviolations = 5;
% percentage spikes missing parameters 
param.maxPercSpikesMissing = 20;
param.computeTimeChunks = 0;
param.deltaTimeChunk = NaN; 
% number of spikes
param.minNumSpikes = 300;
% waveform parameters
param.maxNPeaks = 2;
param.maxNTroughs = 1;
param.somatic = 1; 
param.minWvDuration = 100; %ms
param.maxWvDuration = 900; %ms
% amplitude parameters
param.rawFolder = [ephysap_path, '/..'];
param.nRawSpikesToExtract = 100; 
param.minAmplitude = 40; 
% recording parametrs
param.ephys_sample_rate = 30000;
param.nChannels = 385;
% distance metric parameters
param.computeDistanceMetrics = 1;
param.nChannelsIsoDist = 4;
param.isoDmin = NaN;
param.lratioMin = NaN;
param.ssMin = NaN; 
% ACG parameters
param.ACGbinSize = 0.001;
param.ACGduration = 1;
% ISI parameters
param.longISI = 2;
% cell classification parameters
param.propISI = 0.1;
param.templateDuration = 400;
param.pss = 40;

%% compute quality metrics 
ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
rerun = 0;
if isempty(qMetricsExist) || rerun 
    [qMetric, unitTypes] = bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx, savePath);
else
    load(fullfile(savePath, 'qMetric.mat'))
    load(fullfile(savePath, 'param.mat'))
    unitType = nan(length(qMetric.percSpikesMissing),1);
    unitType(qMetric.nPeaks > param.maxNPeaks | qMetric.nTroughs > param.maxNTroughs |qMetric.somatic ~= param.somatic) = 0; %NOISE OR NON-SOMATIC
    unitType(qMetric.percSpikesMissing <= param.maxPercSpikesMissing & qMetric.nSpikes > param.minNumSpikes & ...
        qMetric.nPeaks <= param.maxNPeaks & qMetric.nTroughs <= param.maxNTroughs & qMetric.Fp <= param.maxRPVviolations & ...
        qMetric.somatic == param.somatic & qMetric.rawAmplitude > param.minAmplitude) = 1;%SINGLE SEXY UNIT
    unitType(isnan(unitType)) = 2;% MULTI UNIT 
end
end
 
% %get memmap
% bc_getRawMemMap;
% 
% ephysData = struct;
% ephysData.spike_times = spikeTimes;
% ephysData.spike_times_timeline = spikeTimes ./ 30000;
% ephysData.spike_templates = spikeTemplates;
% ephysData.templates = templateWaveforms;
% ephysData.template_amplitudes = templateAmplitudes;
% ephysData.channel_positions = readNPY([ephysPath filesep 'channel_positions.npy']);
% ephysData.ephys_sample_rate = 30000;
% ephysData.waveform_t = 1e3*((0:size(templateWaveforms, 2) - 1) / 30000);
% ephysParams = struct;
% plotRaw = 0;
% probeLocation=[];
% %keep ap_data ephysData qMetric param probeLocation unitType plotRaw animal experiments iDay
% 
% bc_unitQualityGUI(ap_data.Data.data,ephysData,qMetric, param, probeLocation, unitType, plotRaw);

%% get histology ? 

%% look at units responses: passive no movement + bhv spont movement (similar xx) 
day = experiments(iDay).day;
experiment = 1;
isSpikeGlx=0;
loadClusters = 0;
AP_load_experimentJF;

task_spike_timeline = spike_times_timeline;
task_stimOn_times = stimOn_times;
task_wheel_move_time = wheel_move_time;
task_wheel_starts = wheel_starts;
task_wheel_types = wheel_types;
task_reward = signals_events.responseTimes(n_trials(1):n_trials(end))';
task_trial_conditions = trial_conditions;

% curr_shank = NaN;
% AP_cellrasterJF({stimOn_times,wheel_move_time,wheel_starts,signals_events.responseTimes(n_trials(1):n_trials(end))'}, ...
% {trial_conditions(:,1),trial_conditions(:,2), wheel_types ...
% trial_conditions(:,3)});

experiment = 2;
AP_load_experimentJF;

passive_spike_timeline = spike_times_timeline;
passive_stimOn_times = stimOn_times;
passive_wheel_move_time = wheel_move_time;
passive_wheel_starts = wheel_starts;
passive_wheel_types = wheel_types;
passive_trial_conditions = trial_conditions;

%psthGUI(task_spike_timeline, task_stimOn_times, task_wheel_move_time, task_wheel_starts, ...
%    task_wheel_types, task_reward, task_trial_conditions, passive_spike_timeline, ...
 %   passive_stimOn_times, passive_wheel_move_time, passive_wheel_starts, ...
 %   passive_wheel_types, passive_trial_conditions)


%get locations (cf script bit andy set on slack)
%plot two raster / cell, all waveforms on side in black, this cell in blue
%(see if it's more wide or more narrow hehe), ACG same (don't do histogram,
%do lines) 
end
%% compare labeling to AP noise manual curation 
% % QQ check correct matching 
% AP_load_experimentJF;
% APnoiseUnits = good_templates == 0;
% removeThese = ~ismember(1:max(spikeTemplates), unique(spikeTemplates));
% APnoiseUnits(removeThese) = [];
% BCbadUnits = goodUnits == 0;
% fracConcordance = sum(BCbadUnits(APnoiseUnits) == 1) / numel(BCbadUnits(APnoiseUnits));
% 
% diffLabeled = find(BCbadUnits(APnoiseUnits) == 1 );
% %plot the waveform of units not indentified as noise by me 
% for iDiffLabeledUnit = 1:length(diffLabeled)
%     figure();
%     minWv = max([-2, -qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)) + 1]);
%     maxWv = min([6-abs(minWv), size(templateWaveforms,3) - qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))]);
%     waveformSelect = abs(maxWv)-6:1:maxWv;
%     yLim = [min(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)))), ...
%         max(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))))];
%     
%     for iSubPlot = 1:6
%         subplot(3,2,iSubPlot)
%         plot(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, ...
%             qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))+waveformSelect(iSubPlot)))
%         xlim([0 82])
%         ylim([yLim(1), yLim(2)])
%         box off; makepretty; 
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%     end
% end

%% look through good/bad cells 

%% wide vs narrow 
