clear all;
close all;
mice = {'JF067', 'JF078', 'JF_AL035', 'JF084', 'JF082'};
iMouse = 1;
mouse = mice{iMouse};
savedirs = {'/home/netshare/zinu/JF067/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF_AL035/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF084/UnitMatch/site1/UnitMatch.mat'}; %'/home/netshare/zinu/JF067/UnitMatch_1_11_enny'

%% ~~get match data~~
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

%% ~~active~~

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 5);


animal = mouse;
protocol = ''; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;


waveforms_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
waveforms_raw_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
acg_long_track = nan(size(theseUnits, 2), max_nRecs, 500);
vis_long_track = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
move_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
rew_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);


UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;

for iUnit = 1:size(theseUnits, 2)
    %figure();
    thisUID = UUIDs(theseUnits(iUnit));
    thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueID == thisUID;
    recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);

    for iRecording = 1:size(recordings_unique, 1)
        %try
        thisRecording = recordings_unique(iRecording);
        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording);

        thisUnit_1idx = thisUnit_0idx + 1;
        thisRawWaveformPath = UniqueIDConversion.Path4UnitNPY_noGoodID{GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording};


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
                keepMe(iExperiment) = contains(block.expDef, 'noGo_stage');
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
        experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
        filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');

        waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
        waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
        try
            JF_load_experiment;
            %AP_cellrasterJF
            spike_templates_0dx_unique = unique(spike_templates_0idx);
            if numel(thisUnit_0idx) > 1
                thisUnit_abs = find(ismember(spike_templates_0dx_unique, thisUnit_0idx));

                waveforms_long_track(iUnit, iRecording, :) = nanmean(waveforms(thisUnit_abs, :));

                waveform_long_raw_tmp = zeros(numel(thisUnit_abs), 82);
                waveform_long_raw_enny_tmp = zeros(numel(thisUnit_abs), 82);
                for iUnit = 1:numel(thisUnit_abs)
                    max_chan = find(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2))));
                    waveform_long_raw_tmp(iUnit, :) = waveforms_raw(thisUnit_abs(iUnit), max_chan, :);

                    raw_enny = dir(thisRawWaveformPath);
                    raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

                    % Detrending
                    raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
                    raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
                    raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
                    [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

                    waveform_long_raw_enny_tmp(iUnit, :) = raw_wv(:, max_chan, 1);

                end
                waveforms_raw_long_track(iUnit, iRecording, :) = nanmean(waveform_long_raw_tmp);
                waveforms_raw_long_track_enny(iUnit, iRecording, :) = nanmean(waveform_long_raw_enny_tmp);

            else
                thisUnit_abs = find(thisUnit_0idx == spike_templates_0dx_unique);
                waveforms_long_track(iUnit, iRecording, :) = waveforms(thisUnit_abs, :);
                max_chan = find(max(max(waveforms_raw(thisUnit_abs, :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs, :, :), 2))));
                waveforms_raw_long_track(iUnit, iRecording, :) = waveforms_raw(thisUnit_abs, max_chan, :);

                raw_enny = dir(thisRawWaveformPath);
                raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

                % Detrending
                raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
                raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
                raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
                [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

                waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:, max_chan, 1);

            end
            % get ACG
            theseSpikeTimes = spike_times_timeline(ismember(spike_templates_0idx, thisUnit_0idx));
            [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
                ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
            ACG = acg(:, 1, 1);
            acg_long_track(iUnit, iRecording, :) = ACG(501:1000);

            % get vis
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
            vis_long_track(iUnit, iRecording, 1:size(curr_psth, 1), :) = curr_psth;

            % get move
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));

            move_times = stimOn_times + stim_to_move;
            align_group_b = align_group_b(~isnan(move_times));
            move_times = move_times(~isnan(move_times));

            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, move_times, align_group_b);
            %if size(curr_psth)
            move_long_track(iUnit, iRecording, 1:size(curr_psth, 1), :) = curr_psth;
            %
            %
            % get rew
            %try
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 3), unique(trial_conditions(:, 3)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_feedback, align_group_b);
            rew_long_track(iUnit, iRecording, :, :) = curr_psth;
            %catch
            %end
        catch
        end

    end
end

close all;

%% vis
for iUnit = 1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(331)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        xlabel('time (ms)')
        ylabel('a.u.')
        hold on;
        title('template waveform')

        subplot(332)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time (ms)')
        ylabel('amplitude (uV)')
        title('mean raw waveform')

        subplot(333)
        plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        xlabel('time (s)')
        ylabel('sp/s')
        title('ACG')

        % subplot(424)
        % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time (s)')
        % ylabel('sp/s)')


        subplot(3, 3, [4:6])
        ss = smoothdata(squeeze(vis_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:500))) ./ nanstd(ss(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from stim onset (s)')
        ylabel('zscore')
        title('-90* azimuth stimulus')

        subplot(3, 3, [7:9])
        ss2 = smoothdata(squeeze(vis_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]) .* 1000;
        ss2 = (ss2 - nanmean(ss2(:, 1:500))) ./ nanstd(ss2(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss2, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from stim onset (s)')
        ylabel('zscore')
        title('0* azimuth stimulus')
        %
        % subplot(427)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')
        % subplot(428)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')

        % % subplot(529)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
        % % subplot(5,2,10)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
    end
    %titles
    for iSubplot = 1:5
        % iSubplot=iSubplot+1
        if iSubplot == 4
            subplot(3, 3, [iSubplot:iSubplot + 2])
        elseif iSubplot == 5
            subplot(3, 3, [iSubplot + 2:iSubplot + 4])
        else
            subplot(3, 3, iSubplot)
        end
        colormap(cmap_cols)
        c = colorbar;
        c.Ticks = [0.05:1 / size(theseRecordings, 1):0.95];
        c.TickLabels = [1:size(theseRecordings, 1)];
        c.Title.String = 'Day #';
    end


    prettify_plot;
end

%% move
for iUnit = 1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(331)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        xlabel('time (ms)')
        ylabel('a.u.')
        hold on;
        title('template waveform')

        subplot(332)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time (ms)')
        ylabel('amplitude (uV)')
        title('mean raw waveform')

        subplot(333)
        plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        xlabel('time (s)')
        ylabel('sp/s')
        title('ACG')

        % subplot(424)
        % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time (s)')
        % ylabel('sp/s)')


        subplot(3, 3, [4:6])
        ss = smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:500))) ./ nanstd(ss(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from move onset (s)')
        ylabel('zscore')
        title('-90* azimuth stimulus')

        subplot(3, 3, [7:9])
        ss2 = smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]) .* 1000;
        ss2 = (ss2 - nanmean(ss2(:, 1:500))) ./ nanstd(ss2(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss2, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from move onset (s)')
        ylabel('zscore')
        title('0* azimuth stimulus')
        %
        % subplot(427)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')
        % subplot(428)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')

        % % subplot(529)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
        % % subplot(5,2,10)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
    end
    %titles
    for iSubplot = 1:5
        % iSubplot=iSubplot+1
        if iSubplot == 4
            subplot(3, 3, [iSubplot:iSubplot + 2])
        elseif iSubplot == 5
            subplot(3, 3, [iSubplot + 2:iSubplot + 4])
        else
            subplot(3, 3, iSubplot)
        end
        colormap(cmap_cols)
        c = colorbar;
        c.Ticks = [0.05:1 / size(theseRecordings, 1):0.95];
        c.TickLabels = [1:size(theseRecordings, 1)];
        c.Title.String = 'Day #';
    end


    prettify_plot;
end

%% rew
for iUnit = 1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(331)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        xlabel('time (ms)')
        ylabel('a.u.')
        hold on;
        title('template waveform')

        subplot(332)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time (ms)')
        ylabel('amplitude (uV)')
        title('mean raw waveform')

        subplot(333)
        plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        xlabel('time (s)')
        ylabel('sp/s')
        title('ACG')

        % subplot(424)
        % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time (s)')
        % ylabel('sp/s)')


        subplot(3, 3, [4:6])
        ss = smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:500))) ./ nanstd(ss(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from rew onset (s)')
        ylabel('zscore')
        title('No rew')

        subplot(3, 3, [7:9])
        ss2 = smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]) .* 1000;
        ss2 = (ss2 - nanmean(ss2(:, 1:500))) ./ nanstd(ss2(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss2, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from rew onset (s)')
        ylabel('zscore')
        title('rew')
        %
        % subplot(427)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')
        % subplot(428)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')

        % % subplot(529)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
        % % subplot(5,2,10)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
    end
    %titles
    for iSubplot = 1:5
        % iSubplot=iSubplot+1
        if iSubplot == 4
            subplot(3, 3, [iSubplot:iSubplot + 2])
        elseif iSubplot == 5
            subplot(3, 3, [iSubplot + 2:iSubplot + 4])
        else
            subplot(3, 3, iSubplot)
        end
        colormap(cmap_cols)
        c = colorbar;
        c.Ticks = [0.05:1 / size(theseRecordings, 1):0.95];
        c.TickLabels = [1:size(theseRecordings, 1)];
        c.Title.String = 'Day #';
    end


    prettify_plot;
end

%% ~passive~

%% get waveform, ACG

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 3);


animal = mouse;
protocol = ''; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;


waveforms_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
waveforms_raw_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
acg_long_track = nan(size(theseUnits, 2), max_nRecs, 500);
vis_long_track = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
%move_long_track_pass = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
%rew_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);


UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;

for iUnit = 1:size(theseUnits, 2)
    %figure();
    thisUID = UUIDs(theseUnits(iUnit));
    thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueID == thisUID;
    recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);

    for iRecording = 1:size(recordings_unique, 1)
        %try
        thisRecording = recordings_unique(iRecording);
        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording);

        thisUnit_1idx = thisUnit_0idx + 1;
        thisRawWaveformPath = UniqueIDConversion.Path4UnitNPY_noGoodID{GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID & ...
            UniqueIDConversion.recsesAll == thisRecording};


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
        experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
        filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');

        waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
        waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
        try
            JF_load_experiment;
            %AP_cellrasterJF
            spike_templates_0dx_unique = unique(spike_templates_0idx);
            if numel(thisUnit_0idx) > 1
                thisUnit_abs = find(ismember(spike_templates_0dx_unique, thisUnit_0idx));

                waveforms_long_track(iUnit, iRecording, :) = nanmean(waveforms(thisUnit_abs, :));

                waveform_long_raw_tmp = zeros(numel(thisUnit_abs), 82);
                waveform_long_raw_enny_tmp = zeros(numel(thisUnit_abs), 82);
                for iUnit = 1:numel(thisUnit_abs)
                    max_chan = find(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2))));
                    waveform_long_raw_tmp(iUnit, :) = waveforms_raw(thisUnit_abs(iUnit), max_chan, :);

                    raw_enny = dir(thisRawWaveformPath);
                    raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

                    % Detrending
                    raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
                    raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
                    raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
                    [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

                    waveform_long_raw_enny_tmp(iUnit, :) = raw_wv(:, max_chan, 1);

                end
                waveforms_raw_long_track(iUnit, thisRecording, :) = nanmean(waveform_long_raw_tmp);
                waveforms_raw_long_track_enny(iUnit, thisRecording, :) = nanmean(waveform_long_raw_enny_tmp);

            else
                thisUnit_abs = find(thisUnit_0idx == spike_templates_0dx_unique);
                waveforms_long_track(iUnit, thisRecording, :) = waveforms(thisUnit_abs, :);
                max_chan = find(max(max(waveforms_raw(thisUnit_abs, :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs, :, :), 2))));
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


        catch
        end

    end
end

bhvData = noGoWorld_behavior_debug({mouse});

bhvData.correct_trials = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 1) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 1);
bhvData.rxn = bhvData.stim_to_moveMean(1:size(RecSesPerUID, 1)+1, 1);
bhvData.protocol(cellfun(@isempty, bhvData.protocol)) = {'0'};

bhvData.correct_trials_nogo = bhvData.noGo(1:size(RecSesPerUID, 1)+1, 3) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 3);
bhvData.correct_trials_go2 = bhvData.goLeft(1:size(RecSesPerUID, 1)+1, 2) ./ bhvData.nTrials(1:size(RecSesPerUID, 1)+1, 2);
bhvData.correct_trials(26) = [];
bhvData.rxn(26) = [];
bhvData.protocol(26) = [];
bhvData.correct_trials_nogo(26) = [];
bhvData.correct_trials_go2(26) = [];
for iS = 1:size(RecSesPerUID, 1) + 1
    stage(iS) = str2num(bhvData.protocol{iS}(end));

end
stage(26) = [];

%% vis
close all;
recordings_unique = unique(UniqueIDConversion.recsesAll);
cmap_cols = crameri('batlow', size(recordings_unique, 1));
figure();
for iUnit = 1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));


    for iRecording = 1:size(recordings_unique, 1)
        %thisRecotding = theseRecording(iRecording);
        subplot(331)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        xlabel('time (ms)')
        ylabel('a.u.')
        hold on;
        title('template waveform')

        subplot(332)
        plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time (ms)')
        ylabel('amplitude (uV)')
        title('mean raw waveform')

        subplot(333)
        plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        xlabel('time (s)')
        ylabel('sp/s')
        title('ACG')

        % subplot(424)
        % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time (s)')
        % ylabel('sp/s)')


        subplot(3, 3, [4:6])
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:500))) ./ nanstd(ss(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from stim onset (s)')
        ylabel('zscore')
        title('-90* azimuth stimulus')

        subplot(3, 3, [7:9])
        ss2 = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 2, :))', 'movmean', [10, 50]) .* 1000;
        ss2 = (ss2 - nanmean(ss2(:, 1:500))) ./ nanstd(ss2(:, 1:500));
        plot([-0.50005:0.001:0.9995], ss2, 'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from stim onset (s)')
        ylabel('zscore')
        title('0* azimuth stimulus')
        %
        % subplot(427)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')
        % subplot(428)
        % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % hold on;
        % xlabel('time from move onset (s)')
        % ylabel('sp/s)')

        % % subplot(529)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
        % % subplot(5,2,10)
        % % plot([-0.50005:0.001:0.9995],smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        % % hold on;
    end
    %titles
    for iSubplot = 1:5
        % iSubplot=iSubplot+1
        if iSubplot == 4
            subplot(3, 3, [iSubplot:iSubplot + 2])
        elseif iSubplot == 5
            subplot(3, 3, [iSubplot + 2:iSubplot + 4])
        else
            subplot(3, 3, iSubplot)
        end
        colormap(cmap_cols(theseRecordings, :))
        c = colorbar;
        c.Ticks = [0.05:1 / size(theseRecordings, 1):0.95];
        c.TickLabels = [theseRecordings];
        c.Title.String = 'Day #';
    end


    prettify_plot;
end

%% plot try one: selected example cells + movement / percentage correct
examples = 1:size(theseUnits, 2);
%recordings_unique = unique(UniqueIDConversion.recsesAll);
%cmap_cols = crameri('batlow', size(recordings_unique, 1));
last_day = find(stage > 3, 1, 'first') + 2;
close all;
%recordings_unique = unique(UniqueIDConversion.recsesAll);
cmap_cols = crameri('batlow', last_day);

figure();
for ii = 1:length(examples)
    iUnit = examples(ii);
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));
    clearvars maxVal minVal
    for iRecording = 1:last_day
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
        maxVal(iRecording) = nanmax(ss);
        minVal(iRecording) = nanmin(ss);
    end

    for iRecording = 1:last_day

        %thisRecotding = theseRecording(iRecording);
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        if length(unique(ss)) > 50
            % subplot(331)
            % plot([0+0.0137:82/3000:(82*82/3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            % xlabel('time (ms)')
            % ylabel('a.u.')
            % hold on;
            % title('template waveform')

            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+last_day);
            hold on;
            plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            %xlabel('time (ms)')
            %ylabel('amplitude (uV)')
            %title('mean raw waveform')
            box off;
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+last_day+1);
            hold on;
            plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            set(gca, 'XScale', 'log')
            box off;
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;
            %xlabel('time (s)')
            %ylabel('sp/s')
            %title('ACG')

            % subplot(424)
            % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
            % hold on;
            % xlabel('time (s)')
            % ylabel('sp/s)')


            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+iRecording);
            hold on;
            box off;
            set(gca, 'XTickLabel', [])
            %set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, 300:1250))', 'movmean', [10, 70]) .* 1000;
            ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
            t = -0.50005:0.001:0.9995;
            plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));
            ylim([nanmin(minVal), nanmax(maxVal)])
            %hold on;
            %xlabel('time from stim onset (s)')
            %ylabel('zscore')
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
    %titles
    % for iSubplot=1:5
    %    % iSubplot=iSubplot+1
    %     if iSubplot == 4
    %         subplot(3,3,[iSubplot:iSubplot+2])
    %     elseif iSubplot==5
    %         subplot(3,3,[iSubplot+2:iSubplot+4])
    %     else
    %         subplot(3,3,iSubplot)
    %     end
    % colormap(cmap_cols(theseRecordings,:))
    % c = colorbar;
    % c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    % c.TickLabels = [theseRecordings];
    % c.Title.String = 'Day #';
    % end


    % prettify_plot;
end

%% plot behavior on top


subplot(length(examples)+2, last_day+2, [1:(last_day) - 1]);
hold on;
plot(1:last_day, [NaN; NaN; bhvData.rxn(1:last_day-2)])

%subplot(length(examples)+3,size(recordings_unique, 1)+2,[(size(recordings_unique, 1))+3:(size(recordings_unique, 1))*2+1]); hold on;
%plot(1:size(recordings_unique, 1), [NaN; NaN; bhvData.correct_trials(1:size(recordings_unique, 1)-2)])

%prettify_plot

%% plot population average for every day on bottom


vis_long_track_pass_pop = nan(max_nRecs, 3, 1500);

recordings_unique = unique([UniqueIDConversion.recsesAll]);

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
    experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
    filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');

    waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
    waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
    try
        JF_load_experiment;
        %AP_cellrasterJF
        spike_templates_0dx_unique = unique(spike_templates_0idx);


        % get vis
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            spike_templates_0dx_unique, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
        vis_long_track_pass_pop(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;


    catch
    end

end

for iRecording = 1:last_day - 1
    clearvars maxVal minVal
    for iRecording = 1:last_day
        ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, 1, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
        maxVal(iRecording) = nanmax(ss);
        minVal(iRecording) = nanmin(ss);
    end
end
for iRecording = 1:last_day - 1
    subplot(length(examples)+2, last_day+2, (length(examples) + 1)*(last_day + 2)+iRecording);
    hold on;
    box off;
    set(gca, 'XTickLabel', [])
    %set(gca, 'YTickLabel', [])
    set(gca, 'LooseInset', get(gca, 'TightInset'))
    box off;

    ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, 1, 300:1250))', 'movmean', [10, 70]) .* 1000;
    ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
    t = -0.50005:0.001:0.9995;
    plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));
    ylim([nanmin(minVal), nanmax(maxVal)])

end

%% plot try 2: overlay
examples = [2, 7, 8, 11, 12, 19, 21];
recordings_unique = unique(UniqueIDConversion.recsesAll);
cmap_cols = crameri('batlow', size(recordings_unique, 1));

figure();
for ii = 1:length(examples)
    iUnit = examples(ii);
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));


    for iRecording = 1:size(recordings_unique, 1)

        %thisRecotding = theseRecording(iRecording);
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, :))', 'movmean', [10, 50]) .* 1000;
        if length(unique(ss)) > 50
            % subplot(331)
            % plot([0+0.0137:82/3000:(82*82/3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            % xlabel('time (ms)')
            % ylabel('a.u.')
            % hold on;
            % title('template waveform')

            subplot(length(examples)+3, 3, (ii + 1)*(1 + 2)+2);
            hold on;
            plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            %xlabel('time (ms)')
            %ylabel('amplitude (uV)')
            %title('mean raw waveform')

            subplot(length(examples)+3, 1+2, (ii + 1)*(1 + 2)+1+2);
            hold on;
            plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            set(gca, 'XScale', 'log')
            %xlabel('time (s)')
            %ylabel('sp/s')
            %title('ACG')

            % subplot(424)
            % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
            % hold on;
            % xlabel('time (s)')
            % ylabel('sp/s)')


            subplot(length(examples)+3, 1+2, (ii + 1)*(1 + 2)+1);
            hold on;
            box off;
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, 1, 300:1350))', 'movmean', [10, 50]) .* 1000;
            ss = (ss - nanmean(ss(:, 200:500))) ./ nanstd(ss(:, 200:500));
            t = -0.50005:0.001:0.9995;
            plot(t(300:1350), ss, 'Color', cmap_cols(iRecording, :));
            %hold on;
            %xlabel('time from stim onset (s)')
            %ylabel('zscore')
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
    %titles
    % for iSubplot=1:5
    %    % iSubplot=iSubplot+1
    %     if iSubplot == 4
    %         subplot(3,3,[iSubplot:iSubplot+2])
    %     elseif iSubplot==5
    %         subplot(3,3,[iSubplot+2:iSubplot+4])
    %     else
    %         subplot(3,3,iSubplot)
    %     end
    % colormap(cmap_cols(theseRecordings,:))
    % c = colorbar;
    % c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    % c.TickLabels = [theseRecordings];
    % c.Title.String = 'Day #';
    % end


    prettify_plot;
end

%

%% load data single cell

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 3);


animal = mouse;
protocol = ''; % (this is the name of the Signals protocol)
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
vis_long_track = nan(size(theseUnits, 2), max_nRecs, 3, 1500);

%move_long_track_pass = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
%rew_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);


UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;
loadClusters = 0;
for iRecording = 1:size(recordings_unique, 1)

    %for iRecording = 1:size(recordings_unique, 1)
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
    experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
    filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');
    JF_load_experiment;
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
            for iUnit = 1:numel(thisUnit_abs)
                max_chan = find(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs(iUnit), :, :), 2))));
                waveform_long_raw_tmp(iUnit, :) = waveforms_raw(thisUnit_abs(iUnit), max_chan, :);

                raw_enny = dir(thisRawWaveformPath);
                raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);

                % Detrending
                raw_wv = permute(raw_wv, [2, 1, 3]); %detrend works over columns
                raw_wv = detrend(raw_wv, 1); % Detrend (linearly) to be on the safe side. OVER TIME!
                raw_wv = permute(raw_wv, [2, 1, 3]); % Put back in order
                [~, max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70, :, :), 3)), [], 1));

                waveform_long_raw_enny_tmp(iUnit, :) = raw_wv(:, max_chan, 1);

            end
            waveforms_raw_long_track(iUnit, thisRecording, :) = nanmean(waveform_long_raw_tmp);
            waveforms_raw_long_track_enny(iUnit, thisRecording, :) = nanmean(waveform_long_raw_enny_tmp);

        else
            thisUnit_abs = find(thisUnit_0idx == spike_templates_0dx_unique);
            waveforms_long_track(iUnit, thisRecording, :) = waveforms(thisUnit_abs, :);
            max_chan = find(max(max(waveforms_raw(thisUnit_abs, :, :), 2)) == max(max(max(waveforms_raw(thisUnit_abs, :, :), 2))));
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


        %catch
        %end

    end

end

%% load data pop
vis_long_track_pass_pop = nan(max_nRecs, 3, 1500);

recordings_unique = unique([UniqueIDConversion.recsesAll]);

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
    experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
    filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');

    waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
    waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
    try
        JF_load_experiment;
        [unitType, qMetric] = bc_qualityMetricsPipeline_JF(animal, day, site, [], 1, '', 0, 0, 1);


        %AP_cellrasterJF
        spike_templates_0dx_unique = unique(spike_templates_0idx);
        theseSpikeTemplates_0idx = spike_templates_0dx_unique(unitType == 1);

        % get vis
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            spike_templates_0dx_unique, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
        vis_long_track_pass_pop(thisRecording, 1:size(curr_psth, 1), :) = curr_psth;


    catch
    end

end

%% plot
imgT = 1;
examples = 1:size(theseUnits, 2);
%recordings_unique = unique(UniqueIDConversion.recsesAll);
%cmap_cols = crameri('batlow', size(recordings_unique, 1));
last_day = find(stage > 3, 1, 'first') + 2;
close all;
%recordings_unique = unique(UniqueIDConversion.recsesAll);
cmap_cols = crameri('batlow', last_day);

figure();
for ii = 1:length(examples)
    iUnit = examples(ii);
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));
    theseRecordings = find(~isnan(waveforms_long_track(iUnit, :, 50)));
    clearvars maxVal minVal
    for iRecording = 1:last_day
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
        ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
        maxVal(iRecording) = nanmax(ss);
        minVal(iRecording) = nanmin(ss);
    end

    for iRecording = 1:last_day

        %thisRecotding = theseRecording(iRecording);
        ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, :))', 'movmean', [10, 50]) .* 1000;
        if length(unique(ss)) > 50
            % subplot(331)
            % plot([0+0.0137:82/3000:(82*82/3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            % xlabel('time (ms)')
            % ylabel('a.u.')
            % hold on;
            % title('template waveform')

            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+last_day);
            hold on;
            plot([0 + 0.0137:82 / 3000:(82 * 82 / 3000)], squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            %xlabel('time (ms)')
            %ylabel('amplitude (uV)')
            %title('mean raw waveform')
            box off;
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+last_day+1);
            hold on;
            plot([0.0005:0.001:0.5], squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
            hold on;
            set(gca, 'XScale', 'log')
            box off;
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;
            %xlabel('time (s)')
            %ylabel('sp/s')
            %title('ACG')

            % subplot(424)
            % plot([0.0005:0.001:0.1],squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
            % hold on;
            % xlabel('time (s)')
            % ylabel('sp/s)')


            subplot(length(examples)+2, last_day+2, (ii)*(last_day + 2)+iRecording);
            hold on;
            box off;
            set(gca, 'XTickLabel', [])
            %set(gca, 'YTickLabel', [])
            set(gca, 'LooseInset', get(gca, 'TightInset'))
            box off;

            ss = smoothdata(squeeze(vis_long_track_pass(iUnit, iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
            ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
            t = -0.50005:0.001:0.9995;
            plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));
            ylim([nanmin(minVal), nanmax(maxVal)])
            %hold on;
            %xlabel('time from stim onset (s)')
            %ylabel('zscore')
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
    %titles
    % for iSubplot=1:5
    %    % iSubplot=iSubplot+1
    %     if iSubplot == 4
    %         subplot(3,3,[iSubplot:iSubplot+2])
    %     elseif iSubplot==5
    %         subplot(3,3,[iSubplot+2:iSubplot+4])
    %     else
    %         subplot(3,3,iSubplot)
    %     end
    % colormap(cmap_cols(theseRecordings,:))
    % c = colorbar;
    % c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    % c.TickLabels = [theseRecordings];
    % c.Title.String = 'Day #';
    % end


    % prettify_plot;
end


subplot(length(examples)+2, last_day+2, [1:(last_day) - 1]);
hold on;
plot(1:last_day, [NaN; NaN; bhvData.rxn(1:last_day-2)])


for iRecording = 1:last_day - 1
    clearvars maxVal minVal

    ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
    ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
    maxVal(iRecording) = nanmax(ss);
    minVal(iRecording) = nanmin(ss);

end
for iRecording = 1:last_day - 1
    subplot(length(examples)+2, last_day+2, (length(examples) + 1)*(last_day + 2)+iRecording);
    hold on;
    box off;
    set(gca, 'XTickLabel', [])
    %set(gca, 'YTickLabel', [])
    set(gca, 'LooseInset', get(gca, 'TightInset'))
    box off;

    ss = smoothdata(squeeze(vis_long_track_pass_pop(iRecording, imgT, 300:1250))', 'movmean', [10, 70]) .* 1000;
    ss = (ss - nanmean(ss(:, 1:200))) ./ nanstd(ss(:, 1:200));
    t = -0.50005:0.001:0.9995;
    plot(t(300:1250), ss, 'Color', cmap_cols(iRecording, :));
    ylim([nanmin(minVal), nanmax(maxVal)])

end