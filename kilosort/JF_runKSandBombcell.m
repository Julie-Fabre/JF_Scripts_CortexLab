
%% preprocess for npx 2.0 probes : run kilosort, extract sync and lfp channels
function JF_runKSandBombcell(animal, date, site, recording, rerunKS, rerunQM, rerunEP)
%% check data
experiment = [];
%kilosort exists?
cl_myPaths;


[ephysAPfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_includingCompressed', site, recording);
if size(ephysAPfile, 2) > 1 && iscell(ephysAPfile)%keep only ap
    ephysAPfile = ephysAPfile{1};
end
rootH = ['/media/julie/Expansion', '/data_temp/'];
rootZ = fileparts(ephysAPfile);
pathToYourConfigFile = [dropboxPath, 'MATLAB/onPaths/Kilosort2/configFiles'];
if contains(rootZ, 'experiment')
    saveFile = [rootZ, filesep, '..', filesep, '..', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
elseif contains(rootZ, 'continuous')
    saveFile = [rootZ, filesep, '../../../', filesep, 'kilosort2', filesep, 'site', num2str(site), filesep];
elseif contains(rootZ, 'recording')
    toot = fileparts(rootZ);
    saveFile = [rootZ, filesep, '..', filesep, '..', filesep, 'kilosort2', filesep,  'site', num2str(site), filesep , 'recording', toot(end), filesep ];
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
kilosortExists = dir(fullfile(saveFile, 'spike_times.npy'));

% sync exsists?
syncExists = dir(fullfile(saveFile, 'sync.mat'));

%% decompress if needed 
saveFileFolder = fullfile('/media/julie/Expansion', animal, date, ['site', num2str(site)]);
fN = dir(ephysAPfile);
decompDataFile = [saveFileFolder, filesep, fN.name(1:end-14), '_bc_decompressed', fN.name(end-13:end-8),'.ap.bin'];

if contains(ephysAPfile, '.cbin') && isempty(dir(decompDataFile)) && (isempty(kilosortExists) || rerunKS || isempty(syncExists))
    % is already decompressed ?
        % if not,. re-decompress
        disp(['decompressing data', ephysAPfile ,  ' to ', saveFileFolder, '...'])
        mkdir(saveFileFolder)
        decompDataFile = bc_extractCbinData(ephysAPfile, [], [], [], saveFileFolder, 0);
        disp('  done decompressing')
    rootZ = fileparts(decompDataFile);
    metaFile = strrep(ephysAPfile, '.cbin', '.meta');
    [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);

else
    rootZ = fileparts(ephysAPfile);
    if contains(ephysAPfile, '.cbin')
        metaFile = strrep(ephysAPfile, '.cbin', '.meta');
    else
        metaFile = strrep(ephysAPfile, '.bin', '.meta');
    end
    [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);

end




%% kilosort
if isempty(kilosortExists) || rerunKS
    % get channel map file
    recompute = 1;
    disp(['getting channel map for ', ephysAPfile ,  '...'])
    chanMapFile = JF_imroToChannelMapLoc(channelMapIMRO,metaFile,recompute);
    disp('  done')
    % kilosort
    master_kilosort;
    %main_kilosort3;
end

% delete catGT file if catgt 
% ephysAPfile
%% extract sync channel

if isempty(syncExists)
    disp(['extracting sync for ', ephysAPfile ,  '...'])

    if size(ephysAPfile, 2) == 2 && iscell(ephysAPfile)%keep only ap
        ephysAPfile = ephysAPfile{1};
    end
    isSpikeGlx = contains(ephysAPfile, '_g'); %spike glx (2.0 probes) or open ephys (3A probes)?
    if isSpikeGlx
        [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
        if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
%            syncFT(ephysAPfile, 385, ephysKSfile);
        end
    else
        AP_preprocess_phase3_newOEJF_onlysync(animal,date, site)
    end
     disp(' done')
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

%% run quality metrics and ephys properties
experiment = 1;
protocol = '';
runQM = 1;
rerunEP = 1;
runEP =1;
plotGUI = 0;
region = ''; 
bc_qualityMetricsPipeline_JF(animal, date, site, recording, experiment, protocol, rerunQM, plotGUI, runQM);
%bc_ephysPropertiesPipeline_JF(animal, date, site, recording, experiment, protocol, rerunEP, runEP, region);

end