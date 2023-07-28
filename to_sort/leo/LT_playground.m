
%% LT_playground
load('/home/netshare/zaru/JF086/2022-11-02/2/2022-11-02_2_JF086_Timeline.mat')

% Get whether stim was flickering
stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
if any(stimScreen_idx)
    stimScreen_flicker = max(Timeline.rawDAQData(:, stimScreen_idx)) - ...
        min(Timeline.rawDAQData(:, stimScreen_idx)) > 2;
end

% Get photodiode flips (compensate for screen flicker)
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
stimScreen_on = Timeline.rawDAQData(:, photodiode_idx) > 0.2;
stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
photodiode_trace = Timeline.rawDAQData(stimScreen_on, photodiode_idx) > photodiode_thresh;

photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
    photodiode_idx), 3);
photodiode_diff_thresh = range(Timeline.rawDAQData(:, photodiode_idx)) * 0.2;
photodiode_diff_t = 20; % time (in ms) to get delayed differential
photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
photodiode_diff_filt = [1, zeros(1, photodiode_diff_samples), -1];
photodiode_trace_diff = abs(conv(photodiode_trace_medfilt, photodiode_diff_filt, 'valid')) > ...
    photodiode_diff_thresh;
photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
    photodiode_trace_diff(2:end)) + photodiode_diff_samples + 1;
photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

figure();
clf;
valu = 100000;
title('Photodiode');
hold on;
plot(photodiode_trace_medfilt(1:valu));
hold on;
scatter(photodiode_flip(find(photodiode_flip <= valu)), ones(size(find(photodiode_flip <= valu), 1), 1))
hold on;
plot(Timeline.rawDAQData(1:valu, photodiode_idx))
hold on;
plot(medfilt1(Timeline.rawDAQData(1:valu, ...
    photodiode_idx), 3))
hold on;
pp = photodiode_flip(find(photodiode_flip <= valu));




% figure out when there were skips
load('/home/netshare/zaru/JF086/2022-11-02/2/2022-11-02_2_JF086_Block.mat')
max(block.events.stimITIsValues)

diff(photodiode_flip_times)
% check this: length(photiodiode_flip_times(2:2:end)) == numberStims 
% match this to stim Ons
numberStims = length(block.events.stimOnTimes);

