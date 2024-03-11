
%% ~~ Figure 5 - example data (B-E)~~
%% Load matches
mouse = 'JF067';
SaveDir = '/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat';
load(SaveDir)

GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
[UniqueIDOpt, idx1, idx2] = unique(UniqueID); %UID options
RecSesOpt = unique(recses); %Recording sessions options
RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
RecSesPerUID = cat(2, RecSesPerUID{:});

%% Get mouse's behavior
bhvData = cl_task_performance({mouse});

bhvData.correct_trials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rxn = bhvData.stim_to_moveMean(1:size(RecSesPerUID, 1)+1, 1);
bhvData.protocol(cellfun(@isempty, bhvData.protocol)) = {'0'};

bhvData.correct_trials_nogo = bhvData.noGo(1:size(RecSesPerUID, 1)+1, 3) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 3);
bhvData.correct_trials_go2 = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 2) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 2);
bhvData.rewardAmount = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) .* bhvData.water_amount(1:size(RecSesPerUID, 1)+1, 1);
bhvData.nRewardedTrials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rewardRate = bhvData.rewardAmount ./ bhvData.duration(1:size(RecSesPerUID, 1)+1)';

if iMouse == 1
    bhvData.correct_trials(26) = [];
    bhvData.rxn(26) = [];
    bhvData.protocol(26) = [];
    bhvData.correct_trials_nogo(26) = [];
    bhvData.correct_trials_go2(26) = [];
    for iS = 1:size(RecSesPerUID, 1) + 1
        stage(iS) = str2num(bhvData.protocol{iS}(end));

    end
    stage(26) = [];
else
    for iS = 1:size(RecSesPerUID, 1) + 1
        stage(iS) = str2num(bhvData.protocol{iS}(end));

    end
end

last_thisDate = find(stage > 3, 1, 'first') + 2;

%% Load single cell PSTHs 

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 3);


animal = mouse;
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
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

waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);


UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;
loadClusters = 0;
for iRecording = 1:last_thisDate

    %for iRecording = 1:size(recordings_unique, 1)
    %try
    thisRecording = iRecording;


    site = 1;
    recording = [];
    thisDate = experiments(thisRecording).thisDate;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);
    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = cl_cortexlab_filename(animal, thisDate, exp, 'block');
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
    filename = cl_cortexlab_filename(animal, thisDate, '', 'ephys', 1, '', '');
    try
        cl_load_experiment;
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
        vis_long_track_pass_depth(iUnit, thisRecording, 1) = template_depths(thisUnit_0idx(1));

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

%% Load population average PSTH
vis_long_track_pass_pop = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop2 = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop_std = nan(max_nRecs, 3, 1500);
vis_long_track_pass_pop_std_n_neuron = nan(max_nRecs);
vis_long_track_pass_pop_whole = nan(max_nRecs, 3, 1500);
recordings_unique = unique([UniqueIDConversion.recsesAll]);
loadClusters = 0;
for iRecording = 1:last_thisDate
    %try
    thisRecording = recordings_unique(iRecording);


    site = 1;
    recording = [];
    thisDate = experiments(thisRecording).thisDate;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);
    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = cl_cortexlab_filename(animal, thisDate, exp, 'block');
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
    experiment = experiments(thisRecording).experiment(1); %n_trials == max(n_trials(logical(keepMe))));
    if length(experiment) > 1
        experiment = experiment(2);
    end
    filename = cl_cortexlab_filename(animal, thisDate, '', 'ephys', 1, '', '');
    try
        cl_load_experiment;
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

        [unitType, qMetric] = bc_qualityMetricsPipeline_JF(animal, thisDate, site, [], 1, '', 0, 0, 1);
    catch
        unitType = ones(length(spike_templates_0dx_unique), 1);
    end
    if contains(block.expDef, 'noGo_stage')
        stimTrials_keep = stim_to_move > 1 | isnan(stim_to_move); % | stim_to_move > 0.8;
    else
        stimTrials_keep = isnan(stim_to_move); % | stim_to_move > 0.8;
    end
    curr_shank = NaN;
    %  cl_cellraster({stimOn_times, stimOn_times}, {trial_conditions(:,1), trial_conditions(:,2)})

    clearvars theseSpikeTemplates_0idx
    if iMouse == 1
        theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1); % & template_depths >= 2100 );
        %theseSpikeTemplates_0idx{1} = spike_templates_0dx_unique(unitType == 1 & template_depths >= 1000 & template_depths <= 2000);
        % theseSpikeTemplates_0idx{2}= spike_templates_0dx_unique(unitType == 1 & template_depths >= 2100 & template_depths <= 2600);
        %theseSpikeTemplates_0idx{3} = spike_templates_0dx_unique(unitType == 1 & template_depths >= 2600 & template_depths <= 3900);
        % theseSpikeTemplates_0idx{4} = spike_templates_0dx_unique(unitType == 1 & template_depths >= 3400 & template_depths <= 3900);
    elseif iMouse == 2
        theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1 & template_depths >= 900 & template_depths <= 1500);
    end

    % get vis
    %for iGroup = 1:size(theseSpikeTemplates_0idx,2)
    [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
    %align_group_b = ones(size(stimOn_times,1),1);
    [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
        theseSpikeTemplates_0idx, raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
    vis_long_track_pass_pop(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;
    %end

    clearvars average_per_neuron
    %for iGroup = 1:size(theseSpikeTemplates_0idx,2)
    for iNeuron = 1:length(theseSpikeTemplates_0idx)

        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            theseSpikeTemplates_0idx(iNeuron), raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
        average_per_neuron(iNeuron, 1:size(curr_psth, 1), :) = curr_psth;
    end
    %  vis_long_track_pass_pop2(thisRecording, 1:size(curr_psth, 1), :, iGroup) = nanmean(squeeze(...
    %     average_per_neuron(find(sum(sum(average_per_neuron(:,:,:, iGroup),2),3)>0),:,:)));
    %end


    % if iMouse == 1
    %     theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1);
    % elseif iMouse == 2
    %     theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1 & template_depths >= 900 & template_depths <= 1500);
    % end

    %     % get vis
    %     [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
    %     %align_group_b = ones(size(stimOn_times,1),1);
    %     [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
    %         spike_templates_0dx_unique, raster_window, psth_bin_size, stimOn_times(stimTrials_keep), align_group_b(stimTrials_keep));
    %     vis_long_track_pass_pop_whole(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;
    %
    %     vis_long_track_pass_pop_std(thisRecording, 2, :) = nanstd(smoothdata(squeeze(average_per_neuron(find(...
    %         sum(sum(average_per_neuron(:,2,:),2),3)>0),2,:)), 2,'movmean', [10,70]));
    %      vis_long_track_pass_pop_std(thisRecording, 1, :) = nanstd(smoothdata(squeeze(average_per_neuron(find(...
    %         sum(sum(average_per_neuron(:,1,:),2),3)>0),2,:)), 2,'movmean', [10,70]));
    %     %nanstd(smoothdata(squeeze(average_per_neuron(find(sum(sum(average_per_neuron(:,:,:),2),3)>0),1,:)),[],1), [10,70]);
    %     vis_long_track_pass_pop_std_n_neuron(iRecording) = length(find(sum(sum(average_per_neuron(:,:,:),2),3)>0));
    %     %vis_long_track_pass_pop_std_n_neuron(iRecording) = length(find(sum(sum(average_per_neuron(:,:,:),2),3)>0));
    % %catch
    % %end

end


%% Plot example cells (E) and their waveforms (D) and ACGs
%close all;
last_thisDate = 11;
imgT = 1;
examples = [12, 14, 19, 28]; %[12, 14, 19, 22, 28]
%last_thisDate = 4;
cmap_cols = crameri('batlow', last_thisDate);

figure();


for ii = 1:length(examples)
    iUnit = examples(ii);
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));
    clearvars maxVal minVal
    for iRecording = 2:last_thisDate
        imgT = 1;
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
        imgT = 2;
        ss2 = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss2 = (ss2 - nanmean(ss2(:, 1:200))) ./ nanmean(ss2(:, 1:200));
        maxVal(iRecording) = nanmax([ss, ss2]);
        minVal(iRecording) = nanmin([ss, ss2]);
    end
    if any(~isnan(maxVal))
        for iRecording = 2:last_thisDate

            %thisRecotding = theseRecording(iRecording);
            ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, :))', 'movmean', [10, 50]) .* 1000;
            if length(unique(ss)) > 50
                % subplot(331)
                % plot([0+0.0137:82/3000:(82*82/3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
                % xlabel('time (ms)')
                % ylabel('a.u.')
                % hold on;
                % title('template waveform')

                subplot(length(examples)+3, last_thisDate+2, (ii + 2 )*(last_thisDate + 2)+last_thisDate);
                hold on;

                plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));

                hold on;
                %xlabel('time (ms)')
                %ylabel('amplitude (uV)')
                %title('mean raw waveform')
                box off;
                if ii==1
                xlabel('time (ms)')
                ylabel('amplitude (uV)')
                end

                set(gca, 'LooseInset', get(gca, 'TightInset'))
                box off;

                subplot(length(examples)+3, last_thisDate+2, (ii + 2 )*(last_thisDate + 2)+last_thisDate+1);
                hold on;
                plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
                hold on;
                %set(gca, 'XScale', 'log')
                box off;
                xlim([0, 0.05])
                %set(gca, 'XTickLabel', [])
                %set(gca, 'YTickLabel', [])
                set(gca, 'LooseInset', get(gca, 'TightInset'))
                box off;
                if ii==1
                xlabel('time (s)')
                ylabel('sp/s')
                end
                %title('ACG')

                % subplot(424)
                % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
                % hold on;
                % xlabel('time (s)')
                % ylabel('sp/s)')


                subplot(length(examples)+3, last_thisDate+2, (ii + 2 )*(last_thisDate + 2)+iRecording-1);
                hold on;
                box off;
                set(gca, 'XTickLabel', [])
                %set(gca, 'YTickLabel', [])
                set(gca, 'LooseInset', get(gca, 'TightInset'))
                box off;
                imgT = 1;
                ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
                ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
                t = -0.50005:0.001:0.9995;
                plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));
                hold on;

                imgT = 2;
                ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
                ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
                t = -0.50005:0.001:0.9995;
                plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :), 'LineStyle', '--');
                hold on;

                %plotshaded(t(300:1250), [ss - squeeze(vis_long_track_pass_std(iUnit, iRecording, imgT, 300:1250))'; ...
                %    ss + squeeze(vis_long_track_pass_std(iUnit, iRecording, imgT, 300:1250))'], cmap_cols(iRecording, :))
                line([0, 0], [nanmin(minVal), nanmax(maxVal)])
                ylim([nanmin(minVal), nanmax(maxVal)])
                %hold on;
                if ii==1 && iRecording == 1
                xlabel('time from stim (s)')
                %ylabel([num2str(ii)]) 
                end
                ylabel(['zscore, cell #', num2str(ii)])
                %if iRecording == 2
                %    legend(['example cell #', num2str(ii)])
                %end

                %title('-90* azimuth stimulus')

                % subplot(3,3,[7:9])
                %  ss2 = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000;
                % ss2 = (ss2 - nanmean(ss2(:,1:500))) ./ nanstd(ss2(:,1:500));
                % plot([-0.50005:0.001:0.9995], ss2, 'Color', cmap_cols(iRecording, :));
                % hold on;
                % xlabel('time from stim onset (s)')
                % ylabel('zscore')
                % title('0* azimuth stimulus')
                % %

            end
        end

    end
end

%% Plot population average (C) 
for iRecording = 2:last_thisDate
    iGroup = 1
    subplot(length(examples)+3, last_thisDate+2, (last_thisDate + 2)*(2 + iGroup - 1)+3+iRecording-4);
    hold on;
    box off;
    set(gca, 'XTickLabel', [])
    %set(gca, 'YTickLabel', [])
    set(gca, 'LooseInset', get(gca, 'TightInset'))
    box off;
    imgT = 1;
    ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
    %ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
    ss_std = (squeeze(vis_long_track_pass_pop_std(iRecording, imgT, 300:1250)) .* 1000 - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200)) ./ ...
        sqrt(vis_long_track_pass_pop_std_n_neuron(iRecording)); %  - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200)) ./sqrt();

    ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
    %ss = (ss - normVal(iRecording)) ./ normVal(iRecording);
    t = -0.50005:0.001:0.9995;
    plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));

    imgT = 2;
    ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
    %ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
    ss_std = (squeeze(vis_long_track_pass_pop_std(iRecording, imgT, 300:1250)) .* 1000 - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200)) ./ ...
        sqrt(vis_long_track_pass_pop_std_n_neuron(iRecording)); %  - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200)) ./sqrt();

    ss = (ss - nanmean(ss(:, 1:200))) ./ nanmean(ss(:, 1:200));
    %ss = (ss - normVal(iRecording)) ./ normVal(iRecording);
    t = -0.50005:0.001:0.9995;
    plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :), 'LineStyle', '--');
    %
    %  line([0, 0], [nanmin(minVal), nanmax(maxVal)])
    ylim([-0.1, 4])
    if iRecording==1
    xlabel('time from stim (s)')
    ylabel('zscore, population')
    end
    plotshaded(t(300:1250), [ss - ss_std'; ...
        ss + ss_std'], cmap_cols(iRecording, :))
    %  line([0, 0], [nanmin(minVal), nanmax(maxVal)])
    %end


    % if iRecording ==2
    %     legend('population average')
    % end


    title(['day ', num2str(iRecording)])


end

%% Plot behavioral data (B)
% reaction time 
subplot(length(examples)+3, last_thisDate+2, 1:last_thisDate-1);
hold on;
plot(1:last_thisDate, [NaN; NaN; bhvData.rxn(1:last_thisDate-2)])
for i = 2:11
    nS(i) = nanstd(bhvData.stim_to_move{i}) ./ sqrt(length(bhvData.stim_to_move{i}));
end
errorbar(1:last_thisDate, [NaN; NaN; bhvData.rxn(1:last_thisDate-2)], nS, nS);
ylabel('reaction time (s)');
xlim([1.5, last_thisDate + 0.5])

% fraction correct trials
subplot(length(examples)+3, last_thisDate+2, last_thisDate+3:last_thisDate*2+1);
hold on;
plot(1:last_thisDate, [NaN; NaN; bhvData.correct_trials(1:last_thisDate-2)])
ylabel('frac. correct trials');


% make plot pretty
prettify_plot('XLimits', 'none', 'YLimits', 'none')

% re-set limits 
subplot(length(examples)+3, last_thisDate+2, 1:last_thisDate-1);
xlim([1.5, last_thisDate + 0.5])
subplot(length(examples)+3, last_thisDate+2, last_thisDate+3:last_thisDate*2+1);
xlim([1.5, last_thisDate + 0.5])


%% Validation with stable functional properties (Analogous to figure 4) (F - K?)

MiceOpt = {'JF067', 'JF078', 'JFAL035', 'JF082', 'JF084'}; % Add all mice you want to analyze. 51 = ventricule mostly, don't include


FromDate = datetime("2024-02-15 09:00:00");
UMFiles = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site2/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site3/UnitMatch.mat', ...
    '/home/netshare/zinu/JF082/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF082/UnitMatch/site2/UnitMatch.mat'};

UMFiles = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat'};
% (UMFiles, whichMetric, groupVector, UseKSLabels, pltDayPairFig, overlearning, overlearningFRs)
summaryFunctionalPlots(UMFiles, 'Corr', '', '', '')

%% Density plots: firing change per day vs ISIcorr and Population resp corr

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 2);


animal = mouse;
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
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
for iRecording = 1:last_thisDate

    %for iRecording = 1:size(recordings_unique, 1)
    %try
    thisRecording = iRecording;


    site = 1;
    recording = [];
    thisDate = experiments(thisRecording).thisDate;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);
    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        [block_filename, block_exists] = cl_cortexlab_filename(animal, thisDate, exp, 'block');
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
    filename = cl_cortexlab_filename(animal, thisDate, '', 'ephys', 1, '', '');
    try
        cl_load_experiment;
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
        vis_long_track_pass_depth(iUnit, thisRecording, 1) = template_depths(thisUnit_0idx(1));

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



