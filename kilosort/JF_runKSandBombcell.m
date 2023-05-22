
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_runKSandBombcell(animal, date, site, recording, rerunKS, rerunQM)

experiment = [];
%% convert raw data from .cbin to .bin and get channel map
myPaths;
figure();

[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_includingCompressed', site, recording);
if size(ephysAPfile, 2) > 1 && iscell(ephysAPfile)%keep only ap
    ephysAPfile = ephysAPfile{1};
end
if contains(ephysAPfile, '.cbin')
    % is already decompressed ?
    saveFileFolder = fullfile('/media/julie/Expansion', animal, date, ['site', num2str(site)]);
    fN = dir(ephysAPfile);
    decompDataFile = [saveFileFolder, filesep, fN.name(1:end-14), '_bc_decompressed', fN.name(end-13:end-8),'.ap.bin'];
    if isempty(dir(decompDataFile))
        % if not,. re-decompress
        
        mkdir(saveFileFolder)
        decompDataFile = bc_extractCbinData(ephysAPfile, [], [], [], saveFileFolder, 0);
    end
    rootZ = fileparts(decompDataFile);
    metaFile = strrep(ephysAPfile, '.cbin', '.meta');
    [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);
else
    rootZ = fileparts(ephysAPfile);
    metaFile = strrep(ephysAPfile, '.bin', '.meta');
    [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);

end

%% get channel map file
recompute = 1;

chanMapFile = JF_imroToChannelMapLoc(channelMapIMRO,metaFile,recompute);

%% kilosort
rootH = [extraHDPath, '/data_temp/'];
 pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
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
%pyKS_saveFile = [toot, filesep, 'pykilosort' filesep, 'site', num2str(site), filesep, 'output', filesep];
kilosortExists = dir(fullfile(saveFile, 'spike_times.npy')) ;
if isempty(kilosortExists) || rerunKS
    master_kilosort;
    %main_kilosort3;
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
%bc_qualityMetricsPipeline_JF(animal, date, site, 1, '', rerunQM, 0, 1)
end