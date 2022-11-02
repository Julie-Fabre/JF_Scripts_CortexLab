
%% params
animal = 'AP084';
site = 1;
experiment = 1;
day = '2020-12-11';

%% get raw data
[ephys_path, protocol_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_ap', site);
raw.n_channels = 384; % number of channels
raw.ephys_datatype = 'int16';
raw.kp = strcat(ephys_path(1:end), filesep, '..');
raw.dat = dir([raw.kp, filesep, '*.dat']);
raw.ephys_ap_filename = [raw.dat(1).folder, filesep, raw.dat(1).name];
raw.ap_dat_dir = dir(raw.ephys_ap_filename);
raw.pull_spikeT = -40:41; % number of points to pull for each waveform
raw.microVoltscaling = 0.19499999284744263; %in structure.oebin for openephys, this never changed so hard-coded here-not loading it in.

raw.dataTypeNBytes = numel(typecast(cast(0, raw.ephys_datatype), 'uint8'));
raw.n_samples = raw.ap_dat_dir.bytes / (raw.n_channels * raw.dataTypeNBytes);
raw.ap_data = memmapfile(raw.ephys_ap_filename, 'Format', {raw.ephys_datatype, [raw.n_channels, raw.n_samples], 'data'});
raw.max_pull_spikes = 200; % number of spikes to pull

%% get ephys

AP_load_experimentJF;

%%

%% plot raw data
cCount = cumsum(repmat(1000, size(1:384, 1), 1), 1) * 1000;
timeSecs = 1;

iChunk = 1;
iChunk = iChunk + 1;
timeChunkStart = 5000 + (iChunk - 1) * (timeSecs * 30000);
timeChunk = iChunk * (timeSecs * 30000);

timeChunkStop = iChunk * (timeSecs * 30000);
theseChans = 1:384;
spacer = 1:1000:1000 * length(theseChans);
% t=timeChunkStart:timeChunkStop;
LinePlotReducer(@plot, t, bsxfun(@plus, double(raw.ap_data.data.data(theseChans, ...
    timeChunkStart:timeChunkStop)), spacer'), 'k');
ylim([min(min(bsxfun(@plus, double(raw.ap_data.data.data(theseChans, ...
    timeChunkStart:timeChunkStop)), spacer'))), max(max(bsxfun(@plus, double(raw.ap_data.data.data(theseChans, ...
    timeChunkStart:timeChunkStop)), spacer')))])

%% overlay detected spikes
theseTimesCenter = spike_times_full;
theseTimesCenter = theseTimesCenter(theseTimesCenter > timeChunkStart/ephys_sample_rate);
theseTimesCenter = theseTimesCenter(theseTimesCenter < timeChunkStop/ephys_sample_rate);
if ~isempty(theseTimesCenter)
    theseTimesFull = reshape(theseTimesCenter.*ephys_sample_rate+raw.pull_spikeT, [size(theseTimesCenter, 1) * size(raw.pull_spikeT, 2), 1]);
end

if ~isempty(theseTimesCenter)
    hold on;
    LinePlotReducer(@plot, double(theseTimesFull), bsxfun(@plus, double(raw.ap_data.data.data(theseChans, ...
         uint32(theseTimesFull))), spacer'), 'r');

end
LinePlotExplorer(gcf);
%overlay the spikes