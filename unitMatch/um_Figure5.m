%% Fig 5 Unit Match 
%% ~~ Loading ~~
clear all;
close all;
mice = {'JF067', 'JF078', 'JF_AL035', 'JF084', 'JF082'};
savedirs = {'/home/netshare/zinu/JF067/UnitMatch_copy_20230922/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF_AL035/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF084/UnitMatch/site1/UnitMatch.mat'}; %'/home/netshare/zinu/JF067/UnitMatch_1_11_enny'

%% load in match data example mouse 
iMouse = 1;
mouse = mice{iMouse};
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

%% load in behavior data 
bhvData = noGoWorld_behavior_debug({mouse});

bhvData.correct_trials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rxn = bhvData.stim_to_moveMean(1:size(RecSesPerUID, 1)+1, 1);
bhvData.protocol(cellfun(@isempty, bhvData.protocol)) = {'0'};

bhvData.correct_trials_nogo = bhvData.noGo(1:size(RecSesPerUID, 1)+1, 3) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 3);
bhvData.correct_trials_go2 = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 2) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 2);

%% load in example ACG, Wv, stim-aligned responses 
Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 3);


animal = mouse;
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;


waveforms_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
waveforms_raw_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
acg_long_track = nan(size(theseUnits, 2), max_nRecs, 500);
vis_long_track_pass = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
vis_long_track_pass_std = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
vis_long_track = nan(size(theseUnits, 2), max_nRecs, 3, 1500);

%move_long_track_pass = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
%rew_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);


UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;
loadClusters = 0;
for iRecording = 1:size(unique([UniqueIDConversion.recsesAll]), 1)

    %for iRecording = 1:size(recordings_unique, 1)
    %try
    thisRecording = iRecording;


    site = 1;
    recording = [];
    day = experiments(thisRecording).day;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);
    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, exp, 'block');
        try
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
        continue; %keepMe(1:size(experiments(thisRecording).experiment, 2)) = 1;
    end
    experiment = experiments(thisRecording).experiment(1);
    if length(experiment) > 1
        experiment = 2;
    end
    filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');
    try
        JF_load_experiment;
    catch
        warning('error')
        continue;
    end
    %AP_cellrasterJF
    spike_templates_0dx_unique = unique(spike_templates_0idx);

    for iUnit = 1:size(theseUnits, 2)
        %figure();
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
        %try
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

        % get vis
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
        vis_long_track_pass(iUnit, thisRecording, 1:size(curr_psth, 1), :) = curr_psth;

        for iT = 1:size(curr_psth, 1)
         
        
             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            thisUnit_0idx, raster_window, psth_bin_size, stimOn_times(align_group_b(1:size(stimOn_times, 1)) == iT), align_group_b(align_group_b(1:size(stimOn_times, 1)) == iT));
        %vis_long_track_pass_std(iUnit, thisRecording, 1:size(curr_psth, 1), :) = nanstd(smoothdata(curr_raster, 'gaussian', [10, 70]));
        [N, T] = size(curr_raster); % N = Number_of_Trials, T = Number_of_Time_Points

        % Calculate proportion of trials with a spike for each time point
        p = sum(curr_raster) / N;

        % Calculate variance for each time point
        variance = p .* (1 - p);

        % Calculate standard error for each time point
        %SE = sqrt(variance / N);
        vis_long_track_pass_std(iUnit, thisRecording, iT, :) = variance;
        end

        %catch
        %end

    end

end
last_day =4;
vis_long_track_pass_pop = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop2 = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop_std = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop_std_n_neuron = nan(max_nRecs);
vis_long_track_pass_pop_whole = nan(max_nRecs, 3, 1500);
recordings_unique = unique([UniqueIDConversion.recsesAll]);
loadClusters= 0;
for iRecording = 1:last_day
    %try
    thisRecording = recordings_unique(iRecording);


    site = 1;
    recording = [];
    day = experiments(thisRecording).day;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);
    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, exp, 'block');
        try
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
        continue; %keepMe(1:size(experiments(thisRecording).experiment, 2)) = 1;
    end
    experiment = experiments(thisRecording).experiment(1);%n_trials == max(n_trials(logical(keepMe))));
    if length(experiment) > 1
        experiment = experiment(2);
    end
    filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');
    try
        JF_load_experiment;
    catch
        %experiment
        %JF_load_experiment;
        warning('error')
        continue;
    end
    %try
    waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
    
    waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
    
        spike_templates_0dx_unique = unique(spike_templates_0idx);
        try
            
        [unitType, qMetric] = bc_qualityMetricsPipeline_JF(animal, day, site, [], 1, '', 0, 0, 1);
        catch
            unitType = ones(length( spike_templates_0dx_unique),1);
        end
        if contains(block.expDef, 'noGo_stage')
            stimTrials_keep = stim_to_move > 1 | isnan(stim_to_move); % | stim_to_move > 0.8;
        else
            stimTrials_keep = isnan(stim_to_move); % | stim_to_move > 0.8;
        end

        %AP_cellrasterJF
        if iMouse == 1
            theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1 & template_depths >= 2000 & template_depths <= 2750);
        elseif iMouse == 2
            theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1 & template_depths >= 900 & template_depths <= 1500);
        end

        % get vis
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        %align_group_b = ones(size(stimOn_times,1),1);
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            spike_templates_0dx_unique, raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
        vis_long_track_pass_pop(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;
        
        clearvars average_per_neuron
        for iNeuron = 1:length(theseSpikeTemplates_0idx)
            
             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            spike_templates_0dx_unique(iNeuron), raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
            average_per_neuron(iNeuron, 1:size(curr_psth, 1), :) = curr_psth;
        end
         vis_long_track_pass_pop2(thisRecording, 1:size(curr_psth, 1), :) = nanmean(average_per_neuron(find(sum(sum(average_per_neuron(:,:,:),2),3)>0),:,:));

        if iMouse == 1
            theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1);
        elseif iMouse == 2
            theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1 & template_depths >= 900 & template_depths <= 1500);
        end

        % get vis
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        %align_group_b = ones(size(stimOn_times,1),1);
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            spike_templates_0dx_unique, raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
        vis_long_track_pass_pop_whole(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;
       
        vis_long_track_pass_pop_std(thisRecording, 2, :) = nanstd(smoothdata(squeeze(average_per_neuron(find(sum(sum(average_per_neuron(:,:,:),2),3)>0),2,:)), 2,'movmean', [10,70]));%nanstd(smoothdata(squeeze(average_per_neuron(find(sum(sum(average_per_neuron(:,:,:),2),3)>0),1,:)),[],1), [10,70]);
        vis_long_track_pass_pop_std_n_neuron(iRecording) = length(find(sum(sum(average_per_neuron(:,:,:),2),3)>0));
        %vis_long_track_pass_pop_std_n_neuron(iRecording) = length(find(sum(sum(average_per_neuron(:,:,:),2),3)>0));
    %catch
    %end

end


%% load in cell correlation 
%% load in # cells for 3 mice 
%% load in histology 

%% ~~ Plotting ~~
figure(); 
%% Fig 5A - task structure + timings ; histology 
subplot(6,3,[1,4,7])

%% Fig 5B - reaction times 
subplot(6,3,2)
hold on;
rt_se = nan(last_day,1);
for iDay = 2:4
    rt_se(iDay) = nanstd(bhvData.stim_to_move{iDay})./sqrt(length(bhvData.stim_to_move{iDay}));
end
errorbar(1:last_day,[NaN; NaN; bhvData.rxn(1:last_day-2)],rt_se,rt_se, 'LineWidth',2); 
ylabel('reaction \newline time (s)');
xlim([1.5, last_day+0.5])
xticks([2,3,4])
xticklabels({'(pre-training)', '1', '2'})
xtickangle(360)

%% Fig 5C - fraction correct trials 
subplot(6,3,5)
hold on;
plot(1:last_day, [NaN; NaN; bhvData.correct_trials(1:last_day-2)])
ylabel('frac. \newline correct trials');
xlim([1.5, last_day+0.5])
xticks([2,3,4])
xticklabels({'(pre-training)', '1', '2'})
xtickangle(360)

%% Fig 5D - Example waveforms and ACG 
subplot(6,3,[10,13,16])

imgType = 2;
examples = [6,19,13,18,12];
cmap_cols = crameri('batlow', last_day);


for ii = 1:length(examples)
    iUnit = examples(ii);
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));
    clearvars maxVal minVal
    for iRecording = 1:last_day
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgType, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
        maxVal(iRecording) = nanmax(ss);
        minVal(iRecording) = nanmin(ss);
    end
if any(~isnan(maxVal))
    for iRecording = 2:last_day

        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgType, :))', 'movmean', [10, 50]) .* 1000;
        if length(unique(ss)) > 50

            subplot(length(examples)+3, last_day+2, (ii+2)*(last_day + 2)+last_day);
            hold on;
            plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            box off;
            xlabel('time (ms)')
            ylabel('amplitude (uV)')
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            subplot(length(examples)+3, last_day+2, (ii+2)*(last_day + 2)+last_day+1);
            hold on;
            plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            xlim([0 0.05])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;
            xlabel('time (s)')
            ylabel('sp/s')



        end
    end

end
end



%% Fig 5E - population + cell stim-aligned averages 
subplot(6,3,[8,11,14,17])

%% Fig 5F - cells correlation with population ? 
subplot(6,3,[3,6,9])

%% Fig 5G - # cells tracked over time for 3 mice 
subplot(6,3,[12,15,18])

