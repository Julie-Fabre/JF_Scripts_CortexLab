%% load data 
animals = {'AP100','AP101','AP104','AP105','AP106'};
passiveProtocol = 'AP_lcrGratingPassive';
bhvProtocol = 'AP_stimWheelRight';

animal = animals{1};
experiments = AP_find_experimentsJF(animal, bhvProtocol, true);
experiments = experiments([experiments.ephys]);

day = experiments(1).day;
experiment = experiments(1).experiment;

ephysPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys');
[spikeTimes, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, usedChannels] = bc_loadEphysData(ephysPath);
ephysap_path = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap');
ephysDirPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_dir');
savePath = fullfile(ephysDirPath, 'qMetrics'); 

%% run qmetrics 
param = struct;
param.plotThis = 0;
% refractory period parameters
param.tauR = 0.0010; %refractory period time (s)
param.tauC = 0.0002; %censored period time (s)
param.maxRPVviolations = 0.2;
% percentage spikes missing parameters 
param.maxPercSpikesMissing = 30;
param.computeTimeChunks = 0;
param.deltaTimeChunk = NaN; 
% number of spikes
param.minNumSpikes = 300;
% waveform parameters
param.maxNPeaks = 2;
param.maxNTroughs = 1;
param.axonal = 0; 
% amplitude parameters
param.rawFolder = [ephysap_path, '/..'];
param.nRawSpikesToExtract = 100; 
param.minAmplitude = 20; 
% recording parametrs
param.ephys_sample_rate = 30000;
param.nChannels = 385;
% distance metric parameters
param.computeDistanceMetrics = 0;
param.nChannelsIsoDist = NaN;
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
[qMetric, goodUnits] = bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
    templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,usedChannels, savePath);

%% compare labeling to AP noise manual curation 
% QQ check correct matching 
AP_load_experimentJF;
APnoiseUnits = good_templates == 0;
removeThese = ~ismember(1:max(spikeTemplates), unique(spikeTemplates));
APnoiseUnits(removeThese) = [];
BCbadUnits = goodUnits == 0;
fracConcordance = sum(BCbadUnits(APnoiseUnits) == 1) / numel(BCbadUnits(APnoiseUnits));

diffLabeled = find(BCbadUnits(APnoiseUnits) == 1 );
%plot the waveform of units not indentified as noise by me 
for iDiffLabeledUnit = 1:length(diffLabeled)
    figure();
    minWv = max([-2, -qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)) + 1]);
    maxWv = min([6-abs(minWv), size(templateWaveforms,3) - qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))]);
    waveformSelect = abs(maxWv)-6:1:maxWv;
    yLim = [min(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)))), ...
        max(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))))];
    
    for iSubPlot = 1:6
        subplot(3,2,iSubPlot)
        plot(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, ...
            qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))+waveformSelect(iSubPlot)))
        xlim([0 82])
        ylim([yLim(1), yLim(2)])
        box off; makepretty; 
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
end

%% look through good/bad cells 

%% wide vs narrow 
