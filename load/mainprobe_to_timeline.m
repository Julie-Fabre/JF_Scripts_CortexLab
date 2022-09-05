function [co] = mainprobe_to_timeline(KS_folder, Timeline, ops, expInfo)
% fundciton takes the spike times from masterprobe and timeline, and
% outputs coefficients to fit timeline
% ops should contain  ops.recording_software='SpikeGLX' ops.ephys_folder
% are ops.ephys_name are just naming tags to access nidaq.bin -- look at
% line 56 to change

%KS_folder=kilosort folder of main probe

% Get flipper flips from timeline
inputNames = {Timeline.hw.inputs.name}';


flipperTrace = Timeline.rawDAQData(:, strcmp(inputNames, 'flipper')) > 2;
flipperFlips = sort([strfind(flipperTrace', [0, 1]), strfind(flipperTrace', [1, 0])])' + 1;
figure();
valu = 5000;
plot(Timeline.rawDAQData(1:valu, strcmp(inputNames, 'flipper')))
hold on;
scatter(flipperFlips(flipperFlips < valu), ones(size(find(flipperFlips < valu), 1), 1))

flipperFlipTimesTimeline = Timeline.rawDAQTimestamps(flipperFlips)';
timelinetime = Timeline.rawDAQTimestamps;
recording_software = ops.recording_software;
% ops must contain ops.ephys_folder if recording software is spikeglx

%% read in some information about spike data from kilosort output
% % Read header information
% headerFID = fopen([KS_folder '\dat_params.txt']);
% headerInfo = textscan(headerFID,'%s %s', 'delimiter',{' = '});
% fclose(headerFID);
% headerInfo = [headerInfo{1}'; headerInfo{2}'];
% header = struct(headerInfo{:});
%
% % Load spike data
% ephysSampleRate = str2double(header.apSampleRate);
% spikeTimes = double(readNPY([KS_folder '\spike_times.npy']))./ephysSampleRate;


if contains(recording_software, 'OpenEphys')
    acqLiveSyncIdx = 2;
    flipperSyncIdx = 4;

    sync = load([KS_folder, '\sync.mat']);
    sync = sync.sync;
    if sum(sync == 0) < 100
        warning('check is sync is saved properly!')
    end
    % Get flipper experiment differences by long delays

    %% this section is to get the times when the timeline sends input to the probes. the recirding is typically longer
    % so when timeline is stopped there is a signal? (Idk how it comes)
    % and I need to find this signal (as the signal to the end of
    % experiemnt so that I can extract the flipper times from the probe
    % that are corresponding to times in timeline.

    flipThresh = 1; % time between flips to define experiment gap (s)
    flipTimes = sync(flipperSyncIdx).timestamps;
    flipperStEnIdx = [[1; find(diff(flipTimes) > flipThresh) + 1], [find(diff(flipTimes) > flipThresh); length(flipTimes)]];
    experimentDurations = diff(flipTimes(flipperStEnIdx), [], 2);
    [~, currExpIdx] = min(abs(experimentDurations-Timeline.rawDAQTimestamps(end))); % as you should start the recroding first
    %     currExpIdx = 4;
    flipperFlipTimesFPGA = flipTimes(flipperStEnIdx(currExpIdx, 1):flipperStEnIdx(currExpIdx, 2));


elseif contains(recording_software, 'SpikeGLX')
    root = ops.ephys_folder;
    try
        nidaqfile = [ops.ephys_folder, '/', sprintf('%s_t0.nidq.bin', ops.ephys_name)];
        meta = ReadMeta_GLX(nidaqfile, root); % read metadata for bin file you are trying to manipulate.
        acqlivesyncidx = 2;
        fippersyncidx = 3;
        d = dir(nidaqfile);
        syncNchans = str2double(meta.nSavedChans);
        samplerate = str2double(meta.niSampRate);
        dat_to_voltage_converter = str2double(meta.niAiRangeMax) / 32768; % taken from Bill Karsh's demoSGLX data. Have to multiply bin values with this number to get the true voltage.

        nSamps = d.bytes / 2 / syncNchans;
        mmf = memmapfile(nidaqfile, 'Format', {'int16', [syncNchans, nSamps], 'x'});

        % flipper at PXIE
        flipper = mmf.Data.x(fippersyncidx, :) * dat_to_voltage_converter;
    catch
        load(sprintf('%s//sync.mat', KS_folder));
        if sum(sync == 0) < 100
            warning('check is sync is saved properly!')
        end

        flipper = sync * (5 / 64); %  this is hardcoded because the ranges are not what we expect
        %samplerate=str2double(recdat.meta.imSampRate); % might want to rewrite
        samplerate = 30000; % that is the IMEC sample rate.
    end
    % construct timestamps for flipper
    timestamps = [0:size(flipper, 2) - 1] / samplerate;
    % binarise flipper;
    flipthres = 2; % to correct - should convert to voltage ...
    binary_flipper = flipper > flipthres;
    %
    flipperFlipsPXIE = sort([strfind(binary_flipper, [0, 1]), strfind(binary_flipper, [1, 0])])' + 1;
    flipperFlipTimesPXIE = timestamps(flipperFlipsPXIE)';
    flipperFlipTimesFPGA = flipperFlipTimesPXIE;


    exps = dir(expInfo);
    exps = exps(~ismember({exps.name}, {'.', '..', 'ephys', '*istology'}));
    if length(exps) > 1 % more than one experiment. thsi assumes you run one timeline for whole session
        disp('more than 1 exp, assuming one single timeline recorded them all')
        for iExp = 1:size(exps, 1)
            timeLineFile = dir([expInfo, filesep, exps(iExp).name, '/*Timeline.mat']);
            timelineAll = load([expInfo, filesep, exps(iExp).name, filesep, timeLineFile.name]);
            inputNames = {timelineAll.Timeline.hw.inputs.name}';

            flipperTraceAll = timelineAll.Timeline.rawDAQData(:, strcmp(inputNames, 'flipper')) > 2;
            flipperFlipsAll = sort([strfind(flipperTraceAll', [0, 1]), strfind(flipperTraceAll', [1, 0])])' + 1;
            flipperFlipTimesTimelineAll{iExp, :} = timelineAll.Timeline.rawDAQTimestamps(flipperFlipsAll)';
        end
        flipperFlipTimesTimeline = flipperFlipTimesTimelineAll{ops.exp, :};
        flipThresh = 1;
        if ~isempty(find(diff(flipperFlipTimesPXIE) > flipThresh))
            
            flipperStEnIdx = [[1; find(diff(flipperFlipTimesPXIE) > flipThresh) + 1], [find(diff(flipperFlipTimesPXIE) > flipThresh); length(flipperFlipTimesPXIE)]];
            experimentDurations = diff(flipperFlipTimesPXIE(flipperStEnIdx), [], 2);
            [~, currExpIdx] = min(abs(experimentDurations-Timeline.rawDAQTimestamps(end)));
            flipperFlipTimesFPGA = flipperFlipTimesPXIE(flipperStEnIdx(currExpIdx, 1):flipperStEnIdx(currExpIdx, 2));
            flipperFlipTimesPXIE = flipperFlipTimesPXIE(flipperStEnIdx(currExpIdx, 1):flipperStEnIdx(currExpIdx, 2));
        end
    end
end

%%

% figure; plot(flipperFlipTimesPXIE);
% hold on; plot(flipperStEnIdx(:,1),flipperFlipTimesPXIE(flipperStEnIdx(:,1)),'r.');

%%
% Check that number of flipper flips in timeline matches ephys and
% apply a bunch of corrections if not
numFlipsDiff = abs(diff([length(flipperFlipTimesFPGA), length(flipperFlipTimesTimeline)]));
if numFlipsDiff > 0 && numFlipsDiff < 20
    disp('WARNING = Flipper flip times different in timeline/ephys');
    if diff([length(flipperFlipTimesFPGA), length(flipperFlipTimesTimeline)]) < 20 && length(flipperFlipTimesFPGA) > 500
        disp('Trying to account for missing flips...');
        while length(flipperFlipTimesTimeline) > length(flipperFlipTimesFPGA)
            compareVect = [flipperFlipTimesFPGA - (flipperFlipTimesFPGA(1)), flipperFlipTimesTimeline(1:length(flipperFlipTimesFPGA)) - flipperFlipTimesTimeline(1)];
            errPoint = find(abs(diff(diff(compareVect, [], 2))) > 0.005, 1);
            flipperFlipTimesTimeline(errPoint+2) = [];
            flipperFlipTimesFPGA(errPoint-2:errPoint+2) = [];
            flipperFlipTimesTimeline(errPoint-2:errPoint+2) = [];
        end
        while length(flipperFlipTimesFPGA) > length(flipperFlipTimesTimeline)
            compareVect = [flipperFlipTimesTimeline - (flipperFlipTimesTimeline(1)), flipperFlipTimesFPGA(1:length(flipperFlipTimesTimeline)) - flipperFlipTimesFPGA(1)];
            errPoint = find(abs(diff(diff(compareVect, [], 2))) > 0.005, 1);
            flipperFlipTimesFPGA(errPoint+2) = [];
            flipperFlipTimesFPGA(errPoint-2:errPoint+2) = [];
            flipperFlipTimesTimeline(errPoint-2:errPoint+2) = [];

            if numel(errPoint) == 0
                flipperFlipTimesFPGA = flipperFlipTimesFPGA(1:length(flipperFlipTimesTimeline));
            end

        end
        compareVect = [flipperFlipTimesFPGA - (flipperFlipTimesFPGA(1)), flipperFlipTimesTimeline - flipperFlipTimesTimeline(1)];


        if isempty(find(abs(diff(diff(compareVect, [], 2))) > 0.005, 1));
            fprintf('Success! \n');
            co = robustfit(flipperFlipTimesFPGA, flipperFlipTimesTimeline);
        end
    end
elseif diff([length(flipperFlipTimesPXIE), length(flipperFlipTimesTimeline)]) == 0
    flipperFlipTimesFPGA = flipperFlipTimesPXIE;
    co = robustfit(flipperFlipTimesFPGA, flipperFlipTimesTimeline);
elseif numFlipsDiff == 0
    co = robustfit(flipperFlipTimesFPGA, flipperFlipTimesTimeline);

end

end
