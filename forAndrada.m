%% load JF data
% this will give you 4 variables you need:
% spike_times_timeline = all spike times for all units, in seconds, aligned
%   to the stimulus presentations
% spike_templates = which spike belongs to which unit (eg if you do
%   spike_times_timeline(spike_templates == 5) you get all spikes for unit
%   5
% stimOn_times = stimOnset time (in seconds) for each trial
% stimIDs = which of the 30 natural images where shown for each trial
%% ~~ DMS and PS dataset info ~~
animals = {'JF020', 'JF024'};
days = [1,4];
sites = [1, 1];
probes = [2, 1];
iDataset = 2; % 1 = dorsomedial striatum dataset, 2 = posterior striatum datatset

protocol = 'JF_natural_images';

animal = animals{iDataset};
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

day = experiments(days(iDataset)).day;
experiment = experiments(days(iDataset)).experiment;
experiment = experiment(1);
site = sites(iDataset);
isSpikeGlx = false;
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true;
loadClusters = false;

ephysPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys',sites(iDataset));
%% ~~ Compute quality metrics ~~
if isempty(dir(fullfile(savePath, 'qMetric.mat')))
%% load ephys data to cpmpute quality metrics 
[spikeTimes, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, usedChannels] = bc_loadEphysData(ephysPath);
ephysap_path = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap',sites(iDataset));
ephysDirPath = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_dir',sites(iDataset));
savePath = fullfile(ephysDirPath, 'qMetrics'); 
mkdir(savePath)
%% quality metric params 
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

%% compute an save the quality metrics 
[qMetric, goodUnits] = bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
    templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,usedChannels, savePath);
else
    load(fullfile(savePath, 'qMetric.mat'))
    load(fullfile(savePath, 'param.mat'))
    goodUnits = qMetric.percSpikesMissing <= param.maxPercSpikesMissing & qMetric.nSpikes > param.minNumSpikes & ...
        qMetric.nPeaks <= param.maxNPeaks & qMetric.nTroughs <= param.maxNTroughs & qMetric.Fp <= param.maxRPVviolations & ...
        qMetric.axonal == param.axonal & qMetric.rawAmplitude > param.minAmplitude;
    
end

%% ~~ load whole experiment ~~
locationKeep = 'CP';% keep only striatal cells as per histology 
AP_load_experimentJF;% keep only 'good' cells as per qualityMetrics in pehys data, align ephys data, get stimOnsets


%% ~~ example: looking at cells ~~ 
curr_shank = NaN;
AP_cellrasterJF({stimOn_times(~isnan(stimOn_times))}, {stimIDs(~isnan(stimOn_times))});