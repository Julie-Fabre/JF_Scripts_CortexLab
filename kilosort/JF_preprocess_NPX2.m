
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_preprocess_NPX2(animal, date, chanMapFile, experiment, site, recording, rerunKS, rerunQM)

%% convert raw data from .cbin to .bin 
myPaths;
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
if size(ephysAPfile, 2) > 1 && iscell(ephysAPfile)%keep only ap
    ephysAPfile = ephysAPfile{1};
end
if contains(ephysAPfile, '.cbin')
    saveFileFolder = fullfile('/media/julie/Elements', animal, date, ['site', num2str(site)]);
    bc_extractCbinData(ephysAPfile, [], [], [], saveFileFolder, []);
else
    rootZ = fileparts(ephysAPfile);

end
%% kilosort


rootH = [extraHDPath, '/data_temp/'];
pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
chanMapFilePath = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles', chanMapFile];
if contains(rootZ, 'experiment')
    saveFile = [rootZ, filesep, '..', filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
elseif contains(rootZ, 'continuous')
    saveFile = [rootZ, filesep, '../../../', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
elseif contains(rootZ, 'recording')
    toot = fileparts(rootZ);
    saveFile = [rootZ, filesep, '..', filesep, '..', filesep, 'kilosort2', filesep, 'recording', toot(end), filesep, 'site', num2str(site), filesep];
else
    
    saveFile = [rootZ, filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
end
if contains(rootZ, 'continuous')
    nChannels = 384;
else
    nChannels = 385;
end
saveDir = dir(saveFile);
if isempty(saveDir)
    mkdir(saveFile)
end
pyKS_saveFile = [toot, filesep, 'pykilosort' filesep, 'site', num2str(site), filesep, 'output', filesep];
kilosortExists = dir(fullfile(saveFile, 'spike_times.npy')) || dir(fullfile(pyKS_saveFile, 'spike_times.npy'));
if isempty(kilosortExists) || rerunKS
    master_kilosort;
end

% delete catGT file if catgt 
% ephysAPfile
%% extract sync channel
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
syncExists = dir(fullfile(saveFile, 'sync.mat'));
if isempty(syncExists)
    
    if size(ephysAPfile, 2) == 2 && iscell(ephysAPfile)%keep only ap
        ephysAPfile = ephysAPfile{1};
    end
    isSpikeGlx = contains(ephysAPfile, '_g'); %spike glx (2.0 probes) or open ephys (3A probes)?
    if isSpikeGlx
        [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
        if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
            syncFT(ephysAPfile, 385, ephysKSfile);
        end
    else
        AP_preprocess_phase3_newOEJF_onlysync(animal,date, site)
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
ephysPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
%ephysPath = strrep(ephysPath, 'kilosort2', 'kilosort2');
ephysap_path = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
ephysDirPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_dir',site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
if isempty(qMetricsExist) || rerunQM 
    disp('getting quality metrics ...')
    [spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, goodChannels] = bc_loadEphysData(ephysPath);

    bc_qualityParamValues; 
  

    bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, goodChannels, savePath);
end
try
    tempFile = dir(fullfile([saveDir(1).folder, 'temp_wh.dat']));
catch
    tempFile = dir(fullfile([saveDir.folder, 'temp_wh.dat']));
end
try
    delete(fullfile(tempFile.folder, tempFile.name))
catch
    disp('did not delete temp file')
end
end