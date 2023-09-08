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

%% enny method
keep MatchTable UniqueIDConversion mouse
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
[UniqueIDOpt, idx1, idx2] = unique(UniqueID); %UID options
RecSesOpt = unique(recses); %Recording sessions options
RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
RecSesPerUID = cat(2, RecSesPerUID{:});

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end-17:end-14)}, 1:size(UniqueIDConversion.Path4UnitNPY,2));

%UniqueIDs_rec = arrayfun(@(X) UniqueID(ismember(RecSesOpt, recses(idx2 == X))), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
%RecSesPerUID = cat(2, RecSesPerUID{:});

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 5);


animal = mouse;
protocol = 'noGo'; % (this is the name of the Signals protocol)
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



%UUID(106)

%UUIDS = [378,379,427,448,939,1167,1424,1454,1551,1660,1728,1743,1748,2050,2131,2150,2190];

UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID,2),1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;

for iUnit = 2%1:size(theseUnits, 2)
    figure();
    thisUID = UUIDs(theseUnits(iUnit));
    thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueID == thisUID;
    recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);
    
    for iRecording = 1:size(recordings_unique, 1)
        thisRecording = recordings_unique(iRecording);
        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID &...
            UniqueIDConversion.recsesAll == thisRecording);
       
        thisUnit_1idx = thisUnit_0idx + 1;
        thisPath = UniqueIDConversion.Path4UnitNPY_noGoodID{GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID &...
            UniqueIDConversion.recsesAll == thisRecording};
    
      
            if numel(thisUnit_0idx) > 1
                for iUnit=1:numel(thisUnit_0idx)

                    
                    raw_enny = dir(thisPath);
                    raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);
                     % Detrending
                    raw_wv = permute(raw_wv,[2,1,3]); %detrend works over columns
                    raw_wv = detrend(raw_wv,1); % Detrend (linearly) to be on the safe side. OVER TIME!
                    raw_wv = permute(raw_wv,[2,1,3]);  % Put back in order
                    [~,max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70,:,:),3)),[],1));
        
                    %max_chan = find(max(max(raw_wv(:,:,1),2)) == max(max(max(raw_wv(:,:,1),2))));
                    %waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:,max_chan,1);
                    waveform_long_raw_enny_tmp(iUnit,:) = raw_wv(:,max_chan,1);
                    
                end
                waveforms_raw_long_track_enny(iUnit, iRecording, :) = nanmean(waveform_long_raw_enny_tmp);
        
            else
                
                raw_enny = dir(thisPath);
                raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);
                 % Detrending
                raw_wv = permute(raw_wv,[2,1,3]); %detrend works over columns
                raw_wv = detrend(raw_wv,1); % Detrend (linearly) to be on the safe side. OVER TIME!
                raw_wv = permute(raw_wv,[2,1,3]);  % Put back in order
                [~,max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70,:,:),3)),[],1));
    
                %max_chan = find(max(max(raw_wv(:,:,1),2)) == max(max(max(raw_wv(:,:,1),2))));
                waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:,max_chan,1);
                
                
            end
hold on;plot(squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :)))
                
            
       

    end
end

for iUnit = 2%1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(521)
        plot(squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(522)
        plot(squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
       
        subplot(523)
        plot(squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        subplot(524)
        plot(squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(525)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(526)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(527)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(528)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(529)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(5,2,10)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
    end
%titles 
    for iSubplot=1:10
    subplot(5,2,iSubplot)
    colormap(cmap_cols)
    c = colorbar;
    c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    c.TickLabels = [1:size(theseRecordings, 1)];
    c.Title.String = 'Day #';
    end


    prettify_plot;
end

%% get behavior
if iMouse == 3 %choieworld
    choiceworld_behavior_JF;
    bhvData = bhv;
    bhvData.correct_trials = sum(bhvData.go_left_trials(1:size(RecSesPerUID, 1), 10:11), 2) ./ sum(bhvData.n_trials_condition(1:size(RecSesPerUID, 1), 10:11), 2);
    bhvData.rxn = nanmean(bhvData.stim_rxn_time(1:size(RecSesPerUID, 1), :), 2);
    stage = repmat(1, [1, size(RecSesPerUID, 1)]);
else
    bhvData = noGoWorld_behavior_debug({mouse});
    if iMouse == 1
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
    else

        bhvData.correct_trials = bhvData.goLeft(1:size(RecSesPerUID, 1), 1) ./ bhvData.nTrials(1:size(RecSesPerUID, 1), 1);
        bhvData.rxn = bhvData.stim_to_moveMean(1:size(RecSesPerUID, 1), 1);
        bhvData.protocol(cellfun(@isempty, bhvData.protocol)) = {'0'};
        for iS = 1:size(RecSesPerUID, 1)
            stage(iS) = str2num(bhvData.protocol{iS}(end));
        end
        %s
    end

end

%% plot
cols = bc_colors(5);

figure()
subplot(6, 1, 1)
title(mouse);
hold on;

plot(1:size(RecSesPerUID, 1), stage, 'Color', cols(3, :))
ylabel('task stage #')

subplot(6, 1, 3)
plot(bhvData.correct_trials, 'Color', cols(1, :));
hold on;
if iMouse == 1
    plot(bhvData.correct_trials_nogo, 'Color', cols(4, :))
    plot(bhvData.correct_trials_go2, 'Color', cols(5, :))
    legend('fraction correct go1', 'fraction correct no go', 'fraction correct go2')
else
    ylabel(['frac. correct', newline, 'trials'])
end


subplot(6, 1, 2)
plot(bhvData.rxn, 'Color', cols(2, :)) %overlay shaded
try
    rt_sem = arrayfun(@(x) nanstd(bhvData.stim_to_move{x}), 1:size(RecSesPerUID, 1));
    hold on;
    plotshaded(1:size(RecSesPerUID, 1), [bhvData.rxn(1:size(RecSesPerUID, 1), 1)' - rt_sem; ...
        bhvData.rxn(1:size(RecSesPerUID, 1), 1)' + rt_sem], cols(:, 2))
catch
end
ylabel(['reaction', newline, 'time (s)'])


subplot(6, 1, 4:6)
h = imagesc(RecSesPerUID(:, sum(RecSesPerUID, 1) > 1)');
colormap(flipud(gray))
xlabel('recording day #')
ylabel('# tracked units')
prettify_plot;

%% waveforms and acg for long-tracked units
%MatchTable = TmpFile.MatchTable;
clearvars num_recs pair_ids pair_recs
for iPair = 1:max(MatchTable.UID1)
    these_indices = find(MatchTable.UID1 == iPair & MatchTable.UID2 == iPair);
    num_recs(iPair) = numel(unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]));
    %recs(iPair)
    % if num_recs(iPair) > 0
    % 
    %     pair_ids(iPair,:) = [MatchTable.ID1(these_indices), MatchTable.ID2(these_indices)];
    %     pair_recs(iPair,:) = [MatchTable.RecSes1(these_indices), MatchTable.RecSes2(these_indices)];
    % end
end

uniqueID_sg = unique([MatchTable.UID1, MatchTable.UID2]);
unitTracked = false(size(uniqueID_sg, 1), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 1);
for iPair = 1:size(uniqueID_sg, 1)
    thisPair = uniqueID_sg(iPair);
    these_indices = MatchTable.UID1 == thisPair & MatchTable.UID2 == thisPair;
    recs = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);
    unitTracked(iPair, recs, :) = true;
end
figure();
imagesc(1-unitTracked)
colormap(gray)
xlabel('recording day #')
ylabel('unit #')
title(mouse)
prettify_plot;

%% 

max_nRecs = max(num_recs);
theseUnits = find(num_recs >= 5);


animal = mouse;
protocol = 'noGo'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;

%RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
%RecSesPerUID = cat(2, RecSesPerUID{:});
%uniqueID_good = UniqueIDConversion.UniqueID(GoodId);
%uniqueID_goodUID = arrayfun(@(X) uniqueID_good(idx2 == X), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
%uniqueID_goodUID = cat(2, uniqueID_goodUID{:});
%oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
%oriID_goodUID = arrayfun(@(X) oriID_good(idx2 == X), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in

%oriID_goodUID = cat(2, oriID_goodUID{:});
waveforms_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
waveforms_raw_long_track = nan(size(theseUnits, 2), max_nRecs, 82);
acg_long_track = nan(size(theseUnits, 2), max_nRecs, 500);
vis_long_track = nan(size(theseUnits, 2), max_nRecs, 3, 1500);
move_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
rew_long_track = nan(size(theseUnits, 2), max_nRecs, 2, 1500);
waveforms_raw_long_track_enny = nan(size(theseUnits, 2), max_nRecs, 82);

for iUnit = 1:size(theseUnits, 2)
  these_indices = find(MatchTable.UID1 == theseUnits(iUnit) & MatchTable.UID2 == theseUnits(iUnit));
  recordings_unique = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);

    for iRecording = 1:size(recordings_unique, 1)
        %iRecording = iRecording+1;
        thisRecording = recordings_unique(iRecording);
        thisUnit_0idx = unique(MatchTable.ID1(MatchTable.RecSes1==thisRecording & MatchTable.UID1 == theseUnits(iUnit)));
MatchTable.MatchProb(MatchTable.RecSes1==thisRecording & MatchTable.UID1 == theseUnits(iUnit))
        thisUnit_1idx = thisUnit_0idx + 1;
    
        try
            % load
            %thisRecording = theseRecordings(iRecording);
            %thisRecording = recordings_unique(iRecording+1);
            site = 1;
            recording = [];
            day = experiments(thisRecording).day;
            n_trials = zeros(size(experiments(thisRecording).experiment,2),1);
            for iExperiment = 1:size(experiments(thisRecording).experiment,2)
                exp = experiments(thisRecording).experiment(iExperiment);
                [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, exp, 'block');
                try
                    load(block_filename)
                    if isfield(block.events, 'stim_idValues')
                        n_trials(iExperiment) = length(block.events.stim_idValues);
                    elseif isfield(block.events, 'stimulusOnTimes')
                        n_trials(iExperiment) = length(block.events.stimulusOnTimes);
                    end
                catch
                    n_trials(iExperiment) = NaN;
                end
            end
            experiment = experiments(thisRecording).experiment(n_trials==max(n_trials));
            filename = AP_cortexlab_filenameJF(animal,day,'','ephys',1,'', '');
            waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
            waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
            JF_load_experiment;
            %AP_cellrasterJF
            spike_templates_0dx_unique = unique(spike_templates_0idx);
            thisUnit_abs = find(thisUnit_0idx==spike_templates_0dx_unique);
            %size(waveforms)
            %size(waveforms_raw)
            % day

            % get waveforms
            waveforms_long_track(iUnit, iRecording, :) = waveforms(thisUnit_abs, :);
            %raw_wv = readNPY(UniqueIDConversion.Path4UnitNPY{theseUnits(iUnit)});
            max_chan = find(max(max(waveforms_raw(thisUnit_abs,:,:),2)) == max(max(max(waveforms_raw(thisUnit_abs,:,:),2))));
            waveforms_raw_long_track(iUnit, iRecording, :) = waveforms_raw(thisUnit_abs,max_chan,:);

            raw_enny = dir([filename, filesep, 'RawWaveforms', filesep, 'Unit' num2str(thisUnit_0idx) '_RawSpikes.npy']);
            raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);
            %[~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
            %spikeMap = readNPY(Path4UnitNPY{uid});
            % Detrending
            raw_wv = permute(raw_wv,[2,1,3]); %detrend works over columns
            raw_wv = detrend(raw_wv,1); % Detrend (linearly) to be on the safe side. OVER TIME!
            raw_wv = permute(raw_wv,[2,1,3]);  % Put back in order
            [~,max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70,:,:),3)),[],1));

            %max_chan = find(max(max(raw_wv(:,:,1),2)) == max(max(max(raw_wv(:,:,1),2))));
            waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:,max_chan,1);
       
  

            % get ACG
            theseSpikeTimes = spike_times_timeline(spike_templates_0idx == thisUnit_0idx);
            [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
                ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
            ACG = acg(:, 1, 1);
            acg_long_track(iUnit, iRecording, :) = ACG(501:1000);

            % get vis
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 1), unique(trial_conditions(:, 1)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b);
            vis_long_track(iUnit, iRecording, 1:size(curr_psth,1), :) = curr_psth;

             % get move
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_move, align_group_b);
            move_long_track(iUnit, iRecording, :, :) = curr_psth;


             % get rew
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 3), unique(trial_conditions(:, 3)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_feedback, align_group_b);
            rew_long_track(iUnit, iRecording, :, :) = curr_psth;
        catch
        end

    end
end
 cmap_cols = crameri('batlow', max_nRecs);
for iUnit = 1:size(theseUnits, 2)
   % thisUnit = Or_UniqueID(theseUnits(iUnit));
    %theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
   
    for iRecording = 1:max_nRecs
        subplot(521)
        plot(squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(522)
        plot(squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
       
        subplot(523)
        plot(squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        subplot(524)
        plot(squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(525)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(526)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(527)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(528)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(529)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(5,2,10)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
    end
%titles 
    for iSubplot=1:10
    subplot(5,2,iSubplot)
    colormap(cmap_cols)
    c = colorbar;
    c.Ticks = [0.05:1/max_nRecs:0.95];
    c.TickLabels = [1:max_nRecs];
    c.Title.String = 'Day #';
    end


    prettify_plot;
end



%% enny method
keep MatchTable UniqueIDConversion mouse
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
[UniqueIDOpt, idx1, idx2] = unique(UniqueID); %UID options
RecSesOpt = unique(recses); %Recording sessions options
RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
RecSesPerUID = cat(2, RecSesPerUID{:});

Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end-17:end-14)}, 1:size(UniqueIDConversion.Path4UnitNPY,2));

%UniqueIDs_rec = arrayfun(@(X) UniqueID(ismember(RecSesOpt, recses(idx2 == X))), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
%RecSesPerUID = cat(2, RecSesPerUID{:});

UUIDs = UniqueID(idx1);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
max_nRecs = max(sum(RecSesPerUID));
theseUnits = find(sum(RecSesPerUID) >= 5);


animal = mouse;
protocol = 'noGo'; % (this is the name of the Signals protocol)
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



%UUID(106)

%UUIDS = [378,379,427,448,939,1167,1424,1454,1551,1660,1728,1743,1748,2050,2131,2150,2190];


for iUnit = 2%1:size(theseUnits, 2)
    thisUID = UUIDs(theseUnits(iUnit));
    thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueID == thisUID;
    recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);
    
    for iRecording = 1:size(recordings_unique, 1)
        thisRecording = recordings_unique(iRecording);
        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueID' == thisUID &...
            UniqueIDConversion.recsesAll == thisRecording);
       
        thisUnit_1idx = thisUnit_0idx + 1;
    
        %try
            site = 1;
            recording = [];
            day = experiments(thisRecording).day;
            n_trials = zeros(size(experiments(thisRecording).experiment,2),1);
            keepMe = zeros(size(experiments(thisRecording).experiment,2),1);
            for iExperiment = 1:size(experiments(thisRecording).experiment,2)
                exp = experiments(thisRecording).experiment(iExperiment);
                [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, exp, 'block');
                try
                    load(block_filename)
                    keepMe(iExperiment) = strcmp(block.expDef, 'noGo');
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
            if sum(keepMe)==0
                keepMe(1:size(experiments(thisRecording).experiment,2)) = 1;
            end
            experiment = experiments(thisRecording).experiment(n_trials==max(n_trials(logical(keepMe))));
            filename = AP_cortexlab_filenameJF(animal,day,'','ephys',1,'', '');
            waveforms_raw = readNPY([filename, filesep, 'templates._bc_rawWaveforms.npy']);
            waveforms_raw_peak = readNPY([filename, filesep, 'templates._bc_rawWaveformPeakChannels.npy']);
            JF_load_experiment;
            %AP_cellrasterJF
            spike_templates_0dx_unique = unique(spike_templates_0idx);
            if numel(thisUnit_0idx) > 1
                thisUnit_abs = find(ismember(spike_templates_0dx_unique,thisUnit_0idx));
                % get waveforms
                waveforms_long_track(iUnit, iRecording, :) = nanmean(waveforms(thisUnit_abs, :));
                %raw_wv = readNPY(UniqueIDConversion.Path4UnitNPY{theseUnits(iUnit)});
                 waveform_long_raw_tmp = zeros(numel(thisUnit_abs),82);
                 waveform_long_raw_enny_tmp = zeros(numel(thisUnit_abs),82);
                for iUnit=1:numel(thisUnit_abs)
                    max_chan = find(max(max(waveforms_raw(thisUnit_abs(iUnit),:,:),2)) == max(max(max(waveforms_raw(thisUnit_abs(iUnit),:,:),2))));
                    %waveforms_raw_long_track(iUnit, iRecording, :) = nanmean(waveforms_raw(thisUnit_abs,max_chan,:));
                    waveform_long_raw_tmp(iUnit,:) = waveforms_raw(thisUnit_abs(iUnit),max_chan,:);
                    
                    raw_enny = dir([filename, filesep, 'RawWaveforms', filesep, 'Unit' num2str(thisUnit_0idx(iUnit)) '_RawSpikes.npy']);
                    raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);
                    %[~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
                    %spikeMap = readNPY(Path4UnitNPY{uid});
                    % Detrending
                    raw_wv = permute(raw_wv,[2,1,3]); %detrend works over columns
                    raw_wv = detrend(raw_wv,1); % Detrend (linearly) to be on the safe side. OVER TIME!
                    raw_wv = permute(raw_wv,[2,1,3]);  % Put back in order
                    [~,max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70,:,:),3)),[],1));
        
                    %max_chan = find(max(max(raw_wv(:,:,1),2)) == max(max(max(raw_wv(:,:,1),2))));
                    %waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:,max_chan,1);
                    waveform_long_raw_enny_tmp(iUnit,:) = raw_wv(:,max_chan,1);
                    
                end
                waveforms_raw_long_track(iUnit, iRecording, :) = nanmean(waveform_long_raw_tmp);
                waveforms_raw_long_track_enny(iUnit, iRecording, :) = nanmean(waveform_long_raw_enny_tmp);
        
            else
                thisUnit_abs = find(thisUnit_0idx==1:max(spike_templates_0dx_unique));
                    % get waveforms
                waveforms_long_track(iUnit, iRecording, :) = waveforms(thisUnit_abs, :);
                %raw_wv = readNPY(UniqueIDConversion.Path4UnitNPY{theseUnits(iUnit)});
                max_chan = find(max(max(waveforms_raw(thisUnit_abs,:,:),2)) == max(max(max(waveforms_raw(thisUnit_abs,:,:),2))));
                waveforms_raw_long_track(iUnit, iRecording, :) = waveforms_raw(thisUnit_abs,max_chan,:);
    
                raw_enny = dir([filename, filesep, 'RawWaveforms', filesep, 'Unit' num2str(thisUnit_0idx) '_RawSpikes.npy']);
                raw_wv = readNPY([raw_enny.folder, filesep, raw_enny.name]);
                %[~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
                %spikeMap = readNPY(Path4UnitNPY{uid});
                % Detrending
                raw_wv = permute(raw_wv,[2,1,3]); %detrend works over columns
                raw_wv = detrend(raw_wv,1); % Detrend (linearly) to be on the safe side. OVER TIME!
                raw_wv = permute(raw_wv,[2,1,3]);  % Put back in order
                [~,max_chan] = nanmax(nanmax(abs(nanmean(raw_wv(35:70,:,:),3)),[],1));
    
                %max_chan = find(max(max(raw_wv(:,:,1),2)) == max(max(max(raw_wv(:,:,1),2))));
                waveforms_raw_long_track_enny(iUnit, iRecording, :) = raw_wv(:,max_chan,1);
                subplot(311); hold on;plot(squeeze(waveforms_long_track(iUnit, iRecording, :)))
                subplot(312); hold on;plot(squeeze(waveforms_raw_long_track(iUnit, iRecording, :)))
                subplot(313); hold on;plot(squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :)))
            end
            %size(waveforms)
            %size(waveforms_raw)
            % day

            
       
  

            % get ACG
            theseSpikeTimes = spike_times_timeline(ismember(spike_templates_0idx,thisUnit_0idx));
            [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
                ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
            ACG = acg(:, 1, 1);
            acg_long_track(iUnit, iRecording, :) = ACG(501:1000);

            % get vis
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 1), unique(trial_conditions(:, 1)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times, align_group_b);
            vis_long_track(iUnit, iRecording, 1:size(curr_psth,1), :) = curr_psth;

             % get move
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_move, align_group_b);
            move_long_track(iUnit, iRecording, :, :) = curr_psth;


             % get rew
             try
            [align_group_a, align_group_b] = ismember(trial_conditions(:, 3), unique(trial_conditions(:, 3)));
            [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx, spike_times_timeline, ...
                thisUnit_0idx, raster_window, psth_bin_size, stimOn_times+stim_to_feedback, align_group_b);
            rew_long_track(iUnit, iRecording, :, :) = curr_psth;
             catch
             end
        %catch
        %end

    end
end

for iUnit = 2%1:size(theseUnits, 2)
    thisUnit = Or_UniqueID(theseUnits(iUnit));
    theseRecordings = find(RecSesPerUID(:, theseUnits(iUnit)));

    figure();
    cmap_cols = crameri('batlow', size(theseRecordings, 1));
    for iRecording = 1:size(theseRecordings, 1)
        subplot(521)
        plot(squeeze(waveforms_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(522)
        plot(squeeze(waveforms_raw_long_track_enny(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
       
        subplot(523)
        plot(squeeze(acg_long_track(iUnit, iRecording, :))', 'Color', cmap_cols(iRecording, :));
        hold on;
        set(gca, 'XScale', 'log')
        subplot(524)
        plot(squeeze(acg_long_track(iUnit, iRecording, 1:100))', 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(525)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(526)
        plot(smoothdata(squeeze(vis_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(527)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(528)
        plot(smoothdata(squeeze(move_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;

        subplot(529)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 1, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
        subplot(5,2,10)
        plot(smoothdata(squeeze(rew_long_track(iUnit, iRecording, 2, :))', 'movmean', [10, 50]).*1000, 'Color', cmap_cols(iRecording, :));
        hold on;
    end
%titles 
    for iSubplot=1:10
    subplot(5,2,iSubplot)
    colormap(cmap_cols)
    c = colorbar;
    c.Ticks = [0.05:1/size(theseRecordings, 1):0.95];
    c.TickLabels = [1:size(theseRecordings, 1)];
    c.Title.String = 'Day #';
    end


    prettify_plot;
end

