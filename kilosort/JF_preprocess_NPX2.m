
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_preprocess_NPX2(animal, date, chanMapFile, experiment, site, recording)

%% kilosort
myPaths;
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
rootZ = fileparts(ephysAPfile);
rootH = tempPath;
pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
chanMapFilePath = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles', chanMapFile];
if contains(rootZ, 'experiment')
    saveFile = [rootZ, '/../../kilosort2/site', num2str(site), filesep];
else
    saveFile = [rootZ, '/../kilosort2/site', num2str(site), filesep];
end
mkdir(saveFile)



master_kilosort;

%% extract sync channel
[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = contains(ephysAPfile, '_g'); %spike glx (2.0 probes) or open ephys (3A probes)?
if isSpikeGlx
    [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
    if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
        syncFT(ephysAPfile, 385, ephysKSfile)
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
end