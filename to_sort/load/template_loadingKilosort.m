% reading kilosort output 
% https://github.com/kwikteam/npy-matlab
% https://djoshea.github.io/neuropixel-utils/kilosort/ 


kilosortPath = '/home/netshare/zinu/JF065/2022-01-18/ephys/kilosort2/site1';

%% getting spikes for templates 
spike_times = readNPY(fullfile(kilosortPath, 'spike_times.npy')); % load spike times in samples. 
% you need to divide by sampling rate (usually 30000 samples/second)  to get spike times in seconds 
spike_times_seconds = double(spike_times) ./ 30000; % you need the fucntion double() to convert from integers to floating point values 
spike_template = readNPY(fullfile(kilosortPath, 'spike_templates.npy')); %load spike templates 
spike_template = spike_template + 1; % 0 indexing to 1-indexing (from python to matlab)

theseSpikes = spike_times(spike_template == 1 ); % get the specific spike times for 1 template 


%% getting template waveforms
template_waveform = readNPY(fullfile(kilosortPath, 'templates.npy')); % load template waveforms : number of templates x time x n channels 
thisWaveform = template_waveform(1 , :, :);

%% template amplitudes 
amplitudes = readNPY(fullfile(kilosortPath, 'amplitudes.npy')); % number of templates x time x n channels 
theseAmplitudes = amplitudes(spike_template == 1 );


%% homework
% 1. get the spike times, waveform and amplitudes for template 108 

% 2. get the spike times for template 218

% 3. plot the amplitudes as a function of spike times for template 108

% 4. get the channel with the maximum amplitude waveform for channel 108
% and plot this waveform 