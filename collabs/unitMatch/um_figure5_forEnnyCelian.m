%% for Enny + Celian

%% Recording info
% using iMouse = 1 for now, but JF078, 84 (2 probes) and 82 (3 probes) are also in the same task + all the recs are processed - so they could be integrated. 
clear all;
close all;
mice = {'JF067', 'JF078', 'JF084', 'JF082'};
iMouse = 1;
mouse = mice{iMouse};
savedirs = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat'}; % only JF067 for now

%% Get mouse behavior 
bhvData = cl_task_performance({mouse});

bhvData.correct_trials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rxn = bhvData.stim_to_moveMean(1:size(RecSesPerUID, 1)+1, 1);
bhvData.protocol(cellfun(@isempty, bhvData.protocol)) = {'0'};

bhvData.correct_trials_nogo = bhvData.noGo(1:size(RecSesPerUID, 1)+1, 3) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 3);
bhvData.correct_trials_go2 = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 2) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 2);
bhvData.rewardAmount = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) .* bhvData.water_amount(1:size(RecSesPerUID, 1)+1, 1);
bhvData.nRewardedTrials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rewardRate = bhvData.rewardAmount ./ bhvData.duration(1:size(RecSesPerUID, 1)+1)';


for iS = 1:size(RecSesPerUID, 1) + 1
    stage(iS) = str2num(bhvData.protocol{iS}(end));
end


last_thisDate = find(stage > 3, 1, 'first') + 2;

%% Get match data
SaveDir = savedirs{iMouse};
%[UniqueID, MatchTable] = AssignUniqueID_POSTUM(SaveDir);
load(SaveDir)

%UniqueIDConversion = TmpFile.UniqueIDConversion;
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
[UniqueIDOpt, idx1, idx2] = unique(UniqueID); %UID options
RecSesOpt = unique(recses); %Recording sessions options
RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
RecSesPerUID = cat(2, RecSesPerUID{:});

%% Get PSTHs for matched cells 

% units to get PSTH for 
Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));
UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 1); % get all cells 
UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;

% experiment info 
animal = mouse;
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
experiments = experiments([experiments.ephys]);
loadClusters = 0; % whether to load phy output or not 

% PSTH and ACG parameters
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;

% pre-allocate array space 
waveforms_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
waveforms_raw_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
acg_long_track = nan(size(theseUnits, 2), max_nRecs, 500);
vis_long_track_pass = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);

for iRecording = 1:last_thisDate

    thisRecording = iRecording;
    site = 1;
    recording = [];
    thisDate = experiments(thisRecording).thisDate;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);

    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = cl_cortexlab_filename(animal, thisDate, exp, 'block');
        try % QQ hacky 
            load(block_filename)
            keepMe(iExperiment) = contains(block.expDef, 'choiceworld');
            if isfield(block.events, 'stim_idValues')
                n_trials(iExperiment) = length(block.events.stim_idValues);
            elseif isfield(block.events, 'stimulusOnTimes')
                n_trials(iExperiment) = length(block.events.stimulusOnTimes);
            end
        catch
            n_trials(iExperiment) = NaN;
            keepMe(iExperiment) = false;
        end
    end

    if sum(keepMe) == 0
        continue; 
    end

    experiment = experiments(thisRecording).experiment(1);
    if length(experiment) > 1
        experiment = 2;
    end

    filename = cl_cortexlab_filename(animal, thisDate, '', 'ephys', 1, '', '');

    try % QQ hacky 
        cl_load_experiment;
    catch
        warning('error')
        continue;
    end

    spike_templates_0dx_unique = unique(spike_templates_0idx);

    for iUnit = 1:size(theseUnits, 2)

        thisUID = UUIDs(theseUnits(iUnit));
        thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueID == thisUID;
        recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);
        if ~ismember(iRecording, recordings_unique)
            continue;
        end

        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording);

        thisUnit_1idx = thisUnit_0idx + 1;
        thisRawWaveformPath = UniqueIDConversion.Path4UnitNPY_noGoodID{GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording};
        waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
        waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);

        if numel(thisUnit_0idx) > 1
            thisUnit_abs = find(ismember(spike_templates_0dx_unique, thisUnit_0idx));

            waveforms_long_track(iUnit, iRecording, :) = nanmean(waveforms(thisUnit_abs, :));

            waveform_long_raw_tmp = zeros(numel(thisUnit_abs), 82);
            waveform_long_raw_enny_tmp = zeros(numel(thisUnit_abs), 82);
            for iiUnit = 1:numel(thisUnit_abs)
                max_chan = find(max(max(waveforms_raw(thisUnit_abs(iiUnit), :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs(iiUnit), :, :), 2))));
                if length(max_chan) > 1
                    max_chan = max_chan(1);
                end
                waveform_long_raw_tmp(iiUnit, :) = waveforms_raw(thisUnit_abs(iiUnit), max_chan, :);

                raw_enny = dir(thisRawWaveformPath);
                raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

                % Detrending
                raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
                raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
                raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
                [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

                waveform_long_raw_enny_tmp(iiUnit, :) = raw_wv(:, max_chan, 1);

            end
            waveforms_raw_long_track(iUnit, thisRecording, :) = nanmean(waveform_long_raw_tmp);
            waveforms_raw_long_track_enny(iUnit, thisRecording, :) = nanmean(waveform_long_raw_enny_tmp);

        else
            thisUnit_abs = find(thisUnit_0idx == spike_templates_0dx_unique);
            waveforms_long_track(iUnit, thisRecording, :) = waveforms(thisUnit_abs, :);
            max_chan = find(max(max(waveforms_raw(thisUnit_abs, :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs, :, :), 2))));
            if length(max_chan) > 1
                max_chan = max_chan(1);
            end
            waveforms_raw_long_track(iUnit, thisRecording, :) = waveforms_raw(thisUnit_abs, max_chan, :);

            raw_enny = dir(thisRawWaveformPath);
            raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

            % Detrending
            raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
            raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
            raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
            [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

            waveforms_raw_long_track_enny(iUnit, thisRecording, :) = raw_wv(:, max_chan, 1);

        end

        % get ACG
        theseSpikeTimes = spike_times_timeline(ismember(spike_templates_0idx, thisUnit_0idx));
        [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
            ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
        ACG = acg(:, 1, 1);
        acg_long_track(iUnit, thisRecording, :) = ACG(501:1000);

        % get visual PSTH
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
        vis_long_track_pass(iUnit, thisRecording, 1:size(curr_psth, 1), :) = curr_psth;

    end

end

vis_long_track_pass(:,:,3,:) = []; %remove 90* stimulus, did not always show it. 
%% save data
% unit numbers 
unit_UniqueIDs = UUIDs(theseUnits);
save('/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/collabs/unitMatch/DATA/all_units_uniqueID.mat', 'unit_UniqueIDs')

% PSTHs : unit_number x recording_number x stimulus_type x time (samples). 
%   stimulus_type = 1 means a controlateral stimulus
%   stimulus type = 2 means a central stimulus. 
save('/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/collabs/unitMatch/DATA/all_units_psth.mat', 'vis_long_track_pass')

% PSTH time 
psth_time = raster_window(1)+psth_bin_size/2:psth_bin_size:raster_window(2)-psth_bin_size/2;
save('/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/collabs/unitMatch/DATA/psth_time.mat', 'psth_time')


%% Functional stuff 
% MiceOpt = {'JF067', 'JF078', 'JFAL035', 'JF082', 'JF084'}; % Add all mice you want to analyze. 51 = ventricule mostly, don't include
% 
% 
% FromDate = datetime("2024-02-15 09:00:00");
% UMFiles = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat', ...
%     '/home/netshare/zinu/JF078/UnitMatch/site1/UnitMatch.mat', ...
%     '/home/netshare/zinu/JF078/UnitMatch/site2/UnitMatch.mat', ...
%     '/home/netshare/zinu/JF078/UnitMatch/site3/UnitMatch.mat', ...
%     '/home/netshare/zinu/JF082/UnitMatch/site1/UnitMatch.mat', ...
%     '/home/netshare/zinu/JF082/UnitMatch/site2/UnitMatch.mat'};
% 
% UMFiles = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat'};
% % (UMFiles, whichMetric, groupVector, UseKSLabels, pltDayPairFig, overlearning, overlearningFRs)
% summaryFunctionalPlots(UMFiles, 'Corr', '', '', '', overlearning, overlearningFRs)
