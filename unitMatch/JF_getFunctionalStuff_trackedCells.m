clear all;
close all;
mice = {'JF067', 'JF078', 'JF_AL035', 'JF084', 'JF082'};
iMouse = 1;
mouse = mice{iMouse};
savedirs = {'/home/netshare/zinu/JF067/UnitMatch/UnitMatch.mat', ...
    '/home/netshare/zinu/JF078/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF_AL035/UnitMatch/site1/UnitMatch.mat', ...
    '/home/netshare/zinu/JF084/UnitMatch/site1/UnitMatch.mat'}; %'/home/netshare/zinu/JF067/UnitMatch_1_11_enny'

%% get match data
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

%% get waveform, ACG

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
                keepMe(iExperiment) = strcmp(block.expDef, 'choice');
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
            keepMe(1:size(experiments(thisRecording).experiment, 2)) = 1;
        end
        experiment = experiments(thisRecording).experiment(n_trials == max(n_trials(logical(keepMe))));
        filename = AP_cortexlab_filenameJF(animal, day, '', 'ephys', 1, '', '');

        waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
        waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
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
            theseSpikeTimes = spike_times_timeline(ismember(spike_templates_0idx,thisUnit_0idx));
            [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
                ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
            ACG = acg(:, 1, 1);
            acg_long_track(iUnit, iRecording, :) = ACG(501:1000);

            % get vis
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times,1)));
            vis_long_track(iUnit, iRecording, 1:size(curr_psth,1), :) = curr_psth;

            %  % get move
            % [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
            % [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            %     thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_move, align_group_b);
            % if size(curr_psth)
            % move_long_track(iUnit, iRecording, :, :) = curr_psth;
            % 
            % 
            %  % get rew
            %  try
            % [align_group_a, align_group_b] = ismember(trial_conditions(:, 3), unique(trial_conditions(:, 3)));
            % [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
            %     thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_feedback, align_group_b);
            % rew_long_track(iUnit, iRecording, :, :) = curr_psth;
            %  catch
            %  end
        %catch
        %end

    end
end

close all;

for iUnit = 1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(331)
        plot([0+0.0137:82/3000:(82*82/3000)], squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        xlabel('time (ms)')
        ylabel('a.u.')
        hold on;
        title('template waveform')

        subplot(332)
        plot([0+0.0137:82/3000:(82*82/3000)],squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
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


        subplot(3,3,[4:6])
        ss = smoothdata(squeeze(vis_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000;
        ss = (ss - nanmean(ss(:,1:500))) ./ nanstd(ss(:,1:500));
        plot([-0.50005:0.001:0.9995],ss,'Color', cmap_cols(iRecording, :));
        hold on;
        xlabel('time from stim onset (s)')
        ylabel('zscore')
        title('-90* azimuth stimulus')

        subplot(3,3,[7:9])
         ss2 = smoothdata(squeeze(vis_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000;
        ss2 = (ss2 - nanmean(ss2(:,1:500))) ./ nanstd(ss2(:,1:500));
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
    for iSubplot=1:5
       % iSubplot=iSubplot+1
        if iSubplot == 4
            subplot(3,3,[iSubplot:iSubplot+2])
        elseif iSubplot==5
            subplot(3,3,[iSubplot+2:iSubplot+4])
        else
            subplot(3,3,iSubplot)
        end
    colormap(cmap_cols)
    c = colorbar;
    c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    c.TickLabels = [1:size(theseRecordings, 1)];
    c.Title.String = 'Day #';
    end


    prettify_plot;
end
