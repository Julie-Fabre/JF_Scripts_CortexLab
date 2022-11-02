
%% look at a chunk of raw data: ap and lfp

%% recording info
animal = 'JF017';
day = '2020-10-22';
site = 2;
ephys_sample_rate = 30000;
n_channels = 384;

ephys_ap_filename = AP_cortexlab_filenameJF(animal, day, [], 'ephys_ap', site, []);
chanPath = AP_cortexlab_filenameJF(animal, day, [], 'ephys', site, []);
[data_path, data_path_exists] = AP_cortexlab_filenameJF(animal, day, [], 'ephys_dir', site);

% (get all recordings within site - assume concat at this point)
lfp_recordings = dir([data_path, filesep, 'experiment*']);
lfp_filename = cellfun(@(x) ...
    [data_path, filesep, x, filesep, 'recording1\continuous\Neuropix-3a-100.1\continuous.dat'], ...
    {lfp_recordings.name}, 'uni', false);

chanMap_0idx = readNPY([chanPath, filesep, 'channel_map.npy']);

%% load ap data
ephys_datatype = 'int16';
ap_dat_dir = dir(ephys_ap_filename);
pull_spikeT = -40:41; %number of points to pull for each spike
microVoltscaling = 0.19499999284744263; %in structure.oebin for openephys
used_channels_idx = 0:max(chanMap_0idx) + 1;
dataTypeNBytes = numel(typecast(cast(0, ephys_datatype), 'uint8')); % determine number of bytes per sample
n_samples = ap_dat_dir.bytes / (n_channels * dataTypeNBytes); % Number of samples per channel
ap_data = memmapfile(ephys_ap_filename, 'Format', {ephys_datatype, [n_channels, n_samples], 'data'});
memMapData = ap_data.data.data;

%% plot ap data
close all;
theseTimes = [50, 51]; %time in s
timeChunkStart = theseTimes(1) * ephys_sample_rate;
timeChunkStop = theseTimes(2) * ephys_sample_rate;
chansToPlot = 50:250;
chansToPlot = chansToPlot';

t = timeChunkStart:timeChunkStop;
cCount = cumsum(repmat(1000, size(chansToPlot, 1), 1), 1);
%close all;
figure();
LinePlotReducer(@plot, t/ephys_sample_rate, double(memMapData(chansToPlot, timeChunkStart:timeChunkStop))+double(cCount), 'k');
LinePlotExplorer(gcf);
xlabel('time(s)')
ylabel('channel')
makeprettyLite;
%% load lfp data
ephys_datatype = 'int16';
lfp_filename=lfp_filename{1};
lfp_dat_dir = dir(lfp_filename);
used_channels_idx = 0:max(chanMap_0idx) + 1;
dataTypeNBytes = numel(typecast(cast(0, ephys_datatype), 'uint8')); % determine number of bytes per sample
n_samples = lfp_dat_dir.bytes / (n_channels * dataTypeNBytes); % Number of samples per channel
lfp_data = memmapfile(lfp_filename, 'Format', {ephys_datatype, [n_channels, n_samples], 'data'});
memMapData = lfp_data.data.data;

%% plot lfp data
close all;
theseTimes = [50, 51]; %time in s
timeChunkStart = theseTimes(1) * ephys_sample_rate;
timeChunkStop = theseTimes(2) * ephys_sample_rate;
chansToPlot = 50:250;
chansToPlot = chansToPlot';

t = timeChunkStart:timeChunkStop;
cCount = cumsum(repmat(1000, size(chansToPlot, 1), 1), 1);
%close all;
figure();
LinePlotReducer(@plot, t/ephys_sample_rate, double(memMapData(chansToPlot, timeChunkStart:timeChunkStop))+double(cCount), 'k');
LinePlotExplorer(gcf);
xlabel('time(s)')
ylabel('channel')
makeprettyLite;