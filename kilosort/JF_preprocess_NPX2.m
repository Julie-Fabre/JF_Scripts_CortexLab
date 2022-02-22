
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_preprocess_NPX2(animal, date, chanMapFile, experiment, site, recording, rerunQM)

%% kilosort

myPaths;
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
if size(ephysAPfile, 2) > 1 %keep only ap
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
            syncFT(ephysAPfile, 385, ephysKSfile);
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
if isempty(qMetricsExist) || rerunQM 
    disp('getting quality metrics ...')
    [spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysPath);

    bc_qualityParamValues; 
  

    bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
end
end