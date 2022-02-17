
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_preprocess_NPX2(animal, date, chanMapFile, experiment, site, recording)

%% kilosort

myPaths;
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
if size(ephysAPfile, 2) == 2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
rootZ = fileparts(ephysAPfile);
rootH = tempPath;
pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
chanMapFilePath = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles', chanMapFile];
if contains(rootZ, 'experiment')
    saveFile = [rootZ, filesep, '..', filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
else
    saveFile = [rootZ, filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
end
saveDir = dir(saveFile);
if isempty(saveDir)
    mkdir(saveFile)
end

kilosortExists = dir(fullfile(saveFile, 'spike_times.npy'));
if isempty(kilosortExists)
    master_kilosort;
end

%% extract sync channel
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
syncExists = dir(fullfile(saveFile, 'sync.mat'));
if isempty(syncExists)
    if size(ephysAPfile, 2) == 2 %keep only ap
        ephysAPfile = ephysAPfile{1};
    end
    isSpikeGlx = contains(ephysAPfile, '_g'); %spike glx (2.0 probes) or open ephys (3A probes)?
    if isSpikeGlx
        [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
        if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
            syncFT(ephysAPfile, 385, ephysKSfile)
        end
    end
end

%% extract lfp
% [filename, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'mainfolder', site, recording);
% mainFolder = filename{1}(1:end - 1);
% lfpDir = dir([mainFolder, filesep, 'lfp', filesep, 'lfp.mat']);
% if isempty(lfpDir) %only compute LFP if not already saved on disk
%     %n_channels = header.n_channels;
%     lfp = JF_get_NPX2_LFP(lfpDir);
%     mkdir([mainFolder, filesep, 'lfp'])
%     save([mainFolder, filesep, 'lfp', filesep, 'lfp.mat'], 'lfp', '-v7.3') %save lfp
% end

%% run quality metrics
ephysPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site);

ephysap_path = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site);
ephysDirPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_dir',site);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
%if isempty(qMetricsExist)
    disp('getting quality metrics ...')
    [spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, usedChannels] = bc_loadEphysData(ephysPath);


    param = struct;
    param.plotThis = 0;
    param.plotGlobal = 1;
    % refractory period parameters
    param.tauR = 0.0020; %refractory period time (s)
    param.tauC = 0.0001; %censored period time (s)
    param.maxRPVviolations = 20;
    % percentage spikes missing parameters
    param.maxPercSpikesMissing = 30;
    param.computeTimeChunks = 0;
    param.deltaTimeChunk = NaN;
    % number of spikes
    param.minNumSpikes = 300;
    % waveform parameters
    param.maxNPeaks = 2;
    param.maxNTroughs = 1;
    param.somatic = 1;
    % amplitude parameters
    param.rawFolder = [ephysap_path, '/..'];
    param.nRawSpikesToExtract = 100;
    param.minAmplitude = 20;
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


    bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, savePath);
%end
end