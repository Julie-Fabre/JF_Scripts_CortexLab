
%% Where to save data - CHANGE THESE PATHS
clear all;
close all;
cl_myPaths;
tmpdatafolder = extraHDPath; %'/media/julie/ExtraHD/data_temp'; % temporary folder for temporary decompression of data

%% Information on mice and recording types - CHANGE THESE MOUSE NAMES AND RECORDING TYPES
MiceOpt = {'JF067', 'JF078', 'JFAL035', 'JF082', 'JF084'}; % Add all mice you want to analyze. 51 = ventricule mostly, don't include
RecordingType(ismember(MiceOpt, {'JF067', 'JF078', 'JF_AL035', 'JF082', 'JF084'})) = {'Chronic'}; % 'Acute' or 'Chronic'
miceDays = [1, 29; 1, 15; 1, 20; 1, 31; 3,19]; % 1:11
sites = [1,1; 1,1; 1,1; 1,2; 1,3];
runMe = 1;
saveJF = 1;
for iMouse = [3,5]%1%1%:size(MiceOpt, 2) 

    %% get all raw ephys and kilosort directories - CHANGE THESE PATHS
    theseSites = sites(iMouse,1):sites(iMouse,2);
    for iSite = 1:size(theseSites,2)
        site = theseSites(iSite);
        experiments = AP_find_experimentsJF(MiceOpt{iMouse}, '', true, site); % find all experiments for this mouse
        experiments = experiments([experiments.ephys]); % keep only experiments that have ephys
    
        ephys_dirs = {experiments.ephys_raw_paths};
        kilosort_dirs = arrayfun(@(x) {[experiments(x).ephys_ks_paths, filesep, 'site', num2str(site)]}, 1:size(experiments, 1));
        thisRecordingType = RecordingType{iMouse};
        mouseName = MiceOpt{iMouse};
    
        [filename, file_exists] = AP_cortexlab_filenameJF(MiceOpt{iMouse}, experiments(end).day, '', 'histo_folder', '', '', '');
        SaveDir = fileparts(filename);
       
        if runMe
            %% subselect two first for testing
            ephys_dirs = ephys_dirs(miceDays(iMouse, 1):miceDays(iMouse, 2));
            kilosort_dirs = kilosort_dirs(miceDays(iMouse, 1):miceDays(iMouse, 2));
            day_str = {experiments.day};
            day_str = day_str(miceDays(iMouse, 1):miceDays(iMouse, 2));
           
            % remove any empty ones 
            nonEmptyCells = ~cellfun(@isempty,ephys_dirs) & ~cellfun(@isempty,kilosort_dirs) &...
                arrayfun(@(x) length(kilosort_dirs{x}) > 6, 1:size(kilosort_dirs,2));
            ephys_dirs = ephys_dirs(nonEmptyCells);
            kilosort_dirs = kilosort_dirs(nonEmptyCells);
            day_str = day_str(nonEmptyCells);
        
            %% run unit match
            set(0, 'DefaultFigureVisible', 'off')
            [PrepareClusInfoparams, UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = um_runUnitMatch(kilosort_dirs, ephys_dirs, SaveDir, tmpdatafolder, thisRecordingType, mouseName);
            set(0, 'DefaultFigureVisible', 'on')
            
            %% get all QMs 
            % [unitType, qMetric] = bc_qualityMetricsPipeline_JF(animal, day, site, recording, experiment_num, protocol, rerunQM, plotGUI, runQM)

            allQmetric = struct;
            for iDay = 1:size(ephys_dirs,2)
                %try
                [unitType, qMetric] = bc_qualityMetricsPipeline_JF(mouseName, day_str{iDay}, site, [], 1, '', 0, 0, 1);
                allQMetric(iDay).qMetric = qMetric; 
            end
            for iDay = 1:size(ephys_dirs,2)
                snr(iDay) = nanmedian(allQMetric(iDay).qMetric.signalToNoiseRatio);
                snr_std(iDay) = nanstd(allQMetric(iDay).qMetric.signalToNoiseRatio);
                ampl(iDay)  = nanmedian(allQMetric(iDay).qMetric.rawAmplitude);
                ampl_std(iDay)  = nanstd(allQMetric(iDay).qMetric.rawAmplitude);
                spD(iDay)  = nanmedian(allQMetric(iDay).qMetric.spatialDecaySlope);
                spD_std(iDay)  = nanstd(allQMetric(iDay).qMetric.spatialDecaySlope);
                
                %catch
                %end
            end
            
            figure();
            subplot(311)
            plot(snr)
            ylabel('SNR')
            title(mouseName)
            

            subplot(312)
            plot(ampl)
            ylabel('waveform amplitude (uV)')

            subplot(313)
            plot(spD)
            ylabel('spatial decay slope')
            xlabel('day #')
            prettify_plot;


            %% check output
            EvaluatingUnitMatch([SaveDir, filesep, 'UnitMatch']);
            ComputeFunctionalScores([SaveDir, filesep, 'UnitMatch'], saveJF)
        
%            DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
%            DrawPairsUnitMatch([SaveDir, filesep, 'UnitMatch'], DrawBlind, saveJF);
        
            % Key presses:
            %   Right arrow: next pair
            %   Left arrow: previous pair
            %   Up arrow: label as match
            %   Down arrow: label as non-match
            %   o: label as I don't know (=uncurated)
            %   m: go to first uncurated pair
            %   p: select pair number
            %   s: SAVE
            recompute = 1;
            %FigureFlick([SaveDir, filesep, 'UnitMatch'], 'julie', recompute, saveJF)
            % your labels are saved in MatchTable.<user>
            
            % move dir 
            movefile([SaveDir, filesep, 'UnitMatch'], [SaveDir, filesep, 'UnitMatch', filesep, 'site', num2str(site)])
        else
            load([SaveDir, filesep, 'UnitMatch', filesep, 'UnitMatch.mat'])
        end
    end
end

% %% attempt 1 get matches across 10 first days
% thisThreshold = 0.5;
% % % find cells with pair-wise matches
% % all_matches = nan(max(MatchTable.RecSes1)-1, 1);
% % matchCount = 1;
% % for iRecording = 1:max(MatchTable.RecSes1)-1
% %     clearvars match_pairID
% %     %theseCandidates = MatchTable.RecSes1==iRecording & MatchTable.MatchProb > 0.5;
% %     sessionMatches = find(MatchTable.RecSes1==iRecording &...
% %         MatchTable.MatchProb > thisThreshold);
% %     match_pairID = [double(MatchTable.ID1(sessionMatches)), double(MatchTable.ID2(sessionMatches)),...
% %         double(MatchTable.UID1(sessionMatches)), double(MatchTable.UID2(sessionMatches)),...
% %         MatchTable.MatchProb(sessionMatches)];
% %     excludeThese = match_pairID(:,1) == ...
% %         match_pairID(:,2);
% %     uniqueMatches = unique(match_pairID(~excludeThese,1)); % get unique matches
% %     for iMatch = 1:size(uniqueMatches, 1)
% %         % find the best match QQ
% %         %best_match = find(match_pairID(:,3) == max(match_pairID(match_pairID(:,1) ==...
% %         %    uniqueMatches(iMatch),3)));
% %         idx = ~excludeThese & match_pairID(:,1) == uniqueMatches(iMatch);
% %         best_match = find(match_pairID(:,3) == match_pairID(:,4));
% %         all_matches(iRecording, matchCount) = MatchTable.ID1(sessionMatches(best_match));
% %         all_matches(iRecording, matchCount) = MatchTable.ID2(sessionMatches(best_match));
% %         matchCount = matchCount + 1;
% %     end
% % end
% 
% % find all UIDs equal, ands remove the others
% % excludeThese = MatchTable.ID1 == MatchTable.ID2;
% % all_matches = table;
% % all_matches.uid1 = MatchTable.UID1(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese) ...
% %     & MatchTable.MatchProb(~excludeThese) > thisThreshold));
% % all_matches.day1 = MatchTable.RecSes1(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese) ...
% %     & MatchTable.MatchProb(~excludeThese) > thisThreshold));
% % all_matches.day2 = MatchTable.RecSes2(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese) ...
% %     & MatchTable.MatchProb(~excludeThese) > thisThreshold));
% % % get unique days
% % matchCounts = histc(all_matches, unique(all_matches));
% % find(matchCounts > 1)
% 
% %% attempt 2
% % pairs = [379,2754,2772,1929,514,746,913,939,968,1676,1780,1826,2138,2159,2341,2451,2555,2749,2773,2866,2871,2922,3029];
% MatchTable = TmpFile.MatchTable;
% clearvars num_recs
% for iPair = 1:max(MatchTable.UID1)
%     these_indices = find(MatchTable.UID1 == iPair & MatchTable.UID2 == iPair);
%     num_recs(iPair) = numel(unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]));
% end
% 
% uniqueID_sg = unique([MatchTable.UID1, MatchTable.UID2]);
% unitTracked = false(size(uniqueID_sg,1), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 1);
% for iPair = 1:size(uniqueID_sg,1)
%     thisPair = uniqueID_sg(iPair);
%     these_indices = MatchTable.UID1 == thisPair & MatchTable.UID2 == thisPair;
%     recs = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);
%     unitTracked(iPair, recs, :) = true;
% end
% figure();
% imagesc(1-unitTracked)
% colormap(gray)
% xlabel('recording day #')
% ylabel('unit #')
% title('JF067')
% prettify_plot;
% 
% 
% 
% % plot average passive response 
% %animal = MiceOpt{iMouse};
% protocol = 'JF_choiceworldStimuli'; % (this is the name of the Signals protocol)
% experiments = AP_find_experimentsJF(animal, protocol, true);
% experiments = experiments([experiments.ephys]);
% raster_window = [-0.5, 1];
% psth_bin_size = 0.001;
% 
% warning off;
% thesePairs = find(num_recs >= 4);
% curr_psth_all = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1500);
% zscore_average_all = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1);
% for iPair = 1:size(thesePairs, 2)
%     thisPair = thesePairs(iPair);
%     these_indices = MatchTable.UID1 == thisPair & MatchTable.UID2 == thisPair;
%     recs = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);
% 
%     for iRec = 1:size(recs, 1)
%         try
%             % load recording
%             thisRecording = recs(iRec);
%             site = 1;
%             recording = [];
%             day = experiments(thisRecording).day;
%             experiment = experiments(thisRecording).experiment(1);
% 
%             JF_load_experiment;
% 
%             % get visually aligned activity for unit(s)
%             theseUnits = unique([MatchTable.ID1(these_indices & MatchTable.RecSes1 == thisRecording); ...
%                 MatchTable.ID2(these_indices & MatchTable.RecSes2 == thisRecording)]);
% 
%             [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
%             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx+1, spike_times_timeline, ...
%                 theseUnits, raster_window, psth_bin_size, stimOn_times, align_group_b);
% 
%             % store visual activity
%             curr_psth_all(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = curr_psth; %psth
%             zscore_average_all(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = ...
%                 (nanmean(curr_psth(:,550:750),2) - nanmean(curr_psth(:,1:150),2)) ./ ...
%                 nanstd(curr_psth(:,1:150), [] , 2); % average zscore 
%         catch
% 
%         end
% 
%     end
% end
% figure();
% imagesc(squeeze(nanmean(curr_psth_all(:,:,2, 600),4)))
% figure();
% subplot(121)
% imagesc(zscore_average_all(:,2:end,1))
% subplot(122)
% imagesc(zscore_average_all(:,2:end,2))
% 
% 
% figure();
% iUnit=12
% cmap_cols = crameri('batlow', 29);
% for iRec=1:size(curr_psth_all,2)
% plot(smoothdata(squeeze(curr_psth_all(iUnit, iRec, 1, :)), 'movmean', [20, 70]), 'Color', cmap_cols(iRec,:)); hold on;
% end
% 
% 
% % plot average task response to stim / reward (movement will be a confound
% % here)
% animal = MiceOpt{iMouse};
% protocol = 'noGo'; % (this is the name of the Signals protocol)
% experiments = AP_find_experimentsJF(animal, protocol, true);
% experiments = experiments([experiments.ephys]);
% raster_window = [-0.5, 1];
% psth_bin_size = 0.001;
% 
% warning on;
% thesePairs = find(num_recs >= 4);
% curr_psth_all_task = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1500);
% zscore_average_all_task = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1);
% curr_psth_all_task_move = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1500);
% zscore_average_all_task_move = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1);
% curr_psth_all_task_rew = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1500);
% zscore_average_all_task_rew = nan(size(thesePairs, 2), size(unique([MatchTable.RecSes1, MatchTable.RecSes2]), 1), 3, 1);
% for iPair = 1:size(thesePairs, 2)
%     thisPair = thesePairs(iPair);
%     these_indices = MatchTable.UID1 == thisPair & MatchTable.UID2 == thisPair;
%     recs = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);
% 
%     for iRec = 1:size(recs, 1)
%         try
%             % load recording
%             thisRecording = recs(iRec);
%             site = 1;
%             recording = [];
%             day = experiments(thisRecording).day;
%             experiment = experiments(thisRecording).experiment(1);
% 
%             JF_load_experiment;
% 
%             % get visually aligned activity for unit(s)
%             theseUnits = unique([MatchTable.ID1(these_indices & MatchTable.RecSes1 == thisRecording); ...
%                 MatchTable.ID2(these_indices & MatchTable.RecSes2 == thisRecording)]);
%             %sp_t = unique(spike)
% 
%             % stimOn
%             [align_group_a, align_group_b] = ismember(trial_conditions(:, 1), unique(trial_conditions(:, 1)));
%             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx+1, spike_times_timeline, ...
%                 theseUnits, raster_window, psth_bin_size, stimOn_times, align_group_b);
% 
%             % store visual activity
%             curr_psth_all_task(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = curr_psth; %psth
%             zscore_average_all_task(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = ...
%                 (nanmean(curr_psth(:,550:700),2) - nanmean(curr_psth(:,1:151),2)) ./ ...
%                 nanstd(curr_psth(:,1:151), [] , 2); % average zscore 
% 
%              % move
%             [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
%             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx+1, spike_times_timeline, ...
%                 theseUnits, raster_window, psth_bin_size, stimOn_times + stim_to_move, align_group_b);
% 
%             % store visual activity
%             curr_psth_all_task_move(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = curr_psth; %psth
%             zscore_average_all_task_move(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = ...
%                 (nanmean(curr_psth(:,550:700),2) - nanmean(curr_psth(:,1:151),2)) ./ ...
%                 nanstd(curr_psth(:,1:151), [] , 2); % average zscore 
% 
%              % reward
%             [align_group_a, align_group_b] = ismember(trial_conditions(:, 3), unique(trial_conditions(:, 3)));
%             [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates_0idx+1, spike_times_timeline, ...
%                 theseUnits, raster_window, psth_bin_size, stimOn_times + stim_to_feedback, align_group_b);
% 
%             % store visual activity
%             curr_psth_all_task_rew(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = curr_psth; %psth
%             zscore_average_all_task_rew(iPair, recs(iRec), 1:size(unique(trial_conditions(:, 2)), 1), :) = ...
%                 (nanmean(curr_psth(:,550:700),2) - nanmean(curr_psth(:,1:151),2)) ./ ...
%                 nanstd(curr_psth(:,1:151), [] , 2); % average zscore 
%         catch
% 
%         end
% 
%     end
% end
% figure();
% subplot(331)
% imagesc(zscore_average_all_task(:,:,1))
% subplot(332)
% imagesc(zscore_average_all_task(:,:,2))
% subplot(333)
% imagesc(zscore_average_all_task(:,:,3))
% subplot(334)
% imagesc(zscore_average_all_task_move(:,:,1))
% subplot(335)
% imagesc(zscore_average_all_task_move(:,:,2))
% subplot(336)
% imagesc(zscore_average_all_task_move(:,:,3))
% subplot(337)
% imagesc(zscore_average_all_task_rew(:,:,1))
% subplot(338)
% imagesc(zscore_average_all_task_rew(:,:,2))
% subplot(339)
% imagesc(zscore_average_all_task_rew(:,:,3))
% % other options: reward response? 
% %% attempt 3
% % tmpfile = dir(fullfile(SaveDir, 'UnitMatch', 'UnitMatch.mat'));
% % 
% % load(fullfile(tmpfile.folder, tmpfile.name));
% % MatchTable = TmpFile.MatchTable; %Extract matchtable
% % UMparam = TmpFile.UMparam; % Extract parameters
% % 
% % % Load AUCS
% % if exist(fullfile(SaveDir, 'UnitMatch', 'AUC.mat'))
% %     AUC = load(fullfile(SaveDir, 'UnitMatch', 'AUC.mat'))';
% %     if iMouse == 1
% %         AUCParams = AUC.AUCStruct.ParamNames;
% %         AUCVals = nan(length(AUCParams), length(MiceOpt));
% %     end
% %     AUCVals(:, iMouse) = AUC.AUCStruct.AUC;
% % end
% % 
% % 
% % % Extract groups
% % 
% % WithinIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
% % MatchIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
% % NonMatchIdx = find((MatchTable.ID1 ~= MatchTable.ID2)); % Not the same unit
% % 
% % % Extract cluster information
% % UniqueIDConversion = TmpFile.UniqueIDConversion;
% % if UMparam.GoodUnitsOnly
% %     GoodId = logical(UniqueIDConversion.GoodID);
% % else
% %     GoodId = true(1, length(UniqueIDConversion.GoodID));
% % end
% % UniqueID = UniqueIDConversion.UniqueID(GoodId);
% % OriID = UniqueIDConversion.OriginalClusID(GoodId);
% % OriIDAll = UniqueIDConversion.OriginalClusID;
% % recses = UniqueIDConversion.recsesAll(GoodId);
% % recsesall = UniqueIDConversion.recsesAll;
% % ndays = length(unique(recses));
% % AllKSDir = UMparam.KSDir; %original KS Dir
% % nclus = length(UniqueID);
% % 
% % tmpUM = cat(2, UMTrackingPerformancePerMouse{:});
% % if UseKSLabels
% %     tmpKS = cat(2, KSTrackingPerformancePerMouse{:});
% % end
% % figure('name', 'Tracking Performance')
% % scatter(1:length(MiceOpt), tmpUM(3, :)./tmpUM(3, :), 20, [0, 0, 0], 'filled')
% % hold on
% % if UseKSLabels
% %     scatter(1:length(MiceOpt), tmpKS(2, :)./tmpKS(3, :), 20, [1, 0, 0], 'filled')
% % end
% % scatter(1:length(MiceOpt), tmpUM(2, :)./tmpUM(3, :), 20, [0, 0, 1], 'filled')
% % 
% % set(gca, 'XTick', 1:length(MiceOpt), 'XTickLabel', MiceOpt, 'XTickLabelRotation', 90)
% % if UseKSLabels
% %     legend('maximum possible', 'Kilosort tracked', 'UnitMatch tracked (concatenated)')
% % else
% %     legend('maximum possible', 'UnitMatch tracked (not concatenated)')
% % end
% % xlim([0.5, length(MiceOpt) + 0.5])
% % ylabel('nUnits')
% % makepretty
% % 
% % TrackingPerformance = nan(3, 0); % Difference between recording number, % Tracked units %maximum possibility
% % TrackingPerformanceKS = nan(3, 0); % Difference between recording number, % Tracked units %maximum possibility
% % 
% % MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
% % for did1 = 1:ndays
% %     for did2 = 1:ndays
% %         if did2 <= did1
% %             continue
% %         end
% %         thesedaysidx = find(ismember(recses, [did1, did2]));
% % 
% %         % can possibly only track these many units:
% %         nMax = min([sum(recses(thesedaysidx) == did1), sum(recses(thesedaysidx) == did2)]);
% %         nMatches = length((thesedaysidx)) - length(unique(UniqueID(thesedaysidx)));
% %         TrackingPerformance = cat(2, TrackingPerformance, [did2 - did1, nMatches, nMax]');
% %         nMatches = length((thesedaysidx)) - length(unique(OriID(thesedaysidx)));
% % 
% %     end
% % end
% % 
% % figure('name', ['TrackingPerformance ', UMparam.SaveDir])
% % subplot(1, 2, 1)
% % scatter(TrackingPerformance(3, :), TrackingPerformance(2, :), 20, [0, 0, 0], 'filled')
% % hold on
% % xlims = get(gca, 'xlim');
% % line([0, max(xlims)], [0, max(xlims)], 'color', [0, 0, 0])
% % xlabel('Number Units Available')
% % ylabel('Number Units Tracked')
% % 
% % subplot(1, 2, 2)
% % scatter(TrackingPerformance(1, :), TrackingPerformance(2, :)./TrackingPerformance(3, :), 20, [0, 0, 0], 'filled')
% % xlabel('\Delta Recordings')
% % ylabel('Proportion tracked cells')

