
%% load data
animals = {'AP100', 'AP101', 'AP104', 'AP105', 'AP106'};
passiveProtocol = 'AP_lcrGratingPassive';
bhvProtocol = 'AP_stimWheelRight';
m2CellPvalstim = [];
m2CellPvalmove = [];
m2CellPSTHmove = [];
m2CellPSTHstim = [];
accCellPvalstim = [];
accCellPvalmove = [];
accCellPSTHmove = [];
accCellPSTHstim = [];
for iAnimal = 1:size(animals, 2)
    animal = animals{iAnimal};
    experimentsBhv = AP_find_experimentsJF(animal, bhvProtocol, true);
    experimentsBhv = experimentsBhv([experimentsBhv.ephys]);
    experimentsPass = AP_find_experimentsJF(animal, passiveProtocol, true);
    experimentsPass = experimentsPass([experimentsPass.ephys]);
    for iDay = 1:size(experimentsBhv, 1)
        
            day = experimentsBhv(iDay).day;
            experimentBhv = experimentsBhv(iDay).experiment;
            experimentPass = experimentsPass(iDay).experiment;
            ephysPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys');
            [spikeTimes, spikeTemplates, ...
                templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysPath);
            ephysap_path = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_ap');
            ephysDirPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_dir');
            savePath = fullfile(ephysDirPath, 'qMetrics');

            %% compute quality metrics
            ephysDirPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_dir');
            qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
            rerun = 1;
            if isempty(qMetricsExist) || rerun
                bc_qualityParamValues;
                param.nChannels = 384;
                [qMetric, unitTypes] = bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
            else
                load(fullfile(savePath, 'qMetric.mat'))
                load(fullfile(savePath, 'param.mat'))
                clearvars unitType;
                bc_qualityParamValues;
                bc_getQualityUnitType;
                bc_plotGlobalQualityMetric;
%                 day = experimentsBhv(iDay).day;
%             experiment = experimentBhv;
%             isSpikeGlx = 0;
%             loadClusters = 0;
%             lfp_channel = 'all';
%             loadLFP = 1;
%             ephys_align = 'cortex';
% 
%             AP_load_experiment;
            
                curr_shank = NaN;
              AP_cellrasterJF({stimOn_times}, ...
              {ones(length(stimOn_times), 1)});
            end


            % %get memmap
            % bc_getRawMemMap;
            %
            % ephysData = struct;
            % ephysData.spike_times = spikeTimes;
            % ephysData.spike_times_timeline = spikeTimes ./ 30000;
            % ephysData.spike_templates = spikeTemplates;
            % ephysData.templates = templateWaveforms;
            % ephysData.template_amplitudes = templateAmplitudes;
            % ephysData.channel_positions = channelPositions;
            % ephysData.ephys_sample_rate = 30000;
            % ephysData.waveform_t = 1e3*((0:size(templateWaveforms, 2) - 1) / 30000);
            % ephysParams = struct;
            % plotRaw = 0;
            % probeLocation=[];
            % %keep ap_data ephysData qMetric param probeLocation unitType plotRaw animal experiments iDay
            %
            % bc_unitQualityGUI(memMapData,ephysData,qMetric, param, probeLocation, unitType, plotRaw);

            %% get histology ?

            %% look at units responses: passive no movement + bhv spont movement (similar xx)
%             day = experimentsBhv(iDay).day;
%             experiment = experimentBhv;
%             isSpikeGlx = 0;
%             loadClusters = 0;
%             lfp_channel = 'all';
%             loadLFP = 1;
%             ephys_align = 'cortex';
% 
%             AP_load_experiment;
% 
%             task_spike_timeline = spike_times_timeline;
%             task_stimOn_times = stimOn_times;
%             task_stim_to_move = stim_to_move;
%             task_stim_to_feedback = stim_to_feedback;
%             task_wheel_starts = wheel_starts(wheel_move_iti_idx); %only ITI moves
%             task_reward = signals_events.responseTimes(n_trials(1):n_trials(end))';
% 
%             %  curr_shank = NaN;
%             %  AP_cellrasterJF({stimOn_times}, ...
%             %  {trial_conditions(:,1)});
% 
%             experiment = experimentPass;
%             clearvars spike_times_timeline stimOn_times wheel_starts trial_conditions
%             AP_load_experiment;
% 
%             passive_spike_timeline = spike_times_timeline;
%             passive_stimOn_times = stimOn_times;
%             passive_wheel_starts = wheel_starts;
%             
%             clearvars spike_times_timeline stimOn_times wheel_starts trial_conditions
%             
%             psthGUI(spike_templates, task_spike_timeline, task_stimOn_times, task_wheel_starts,...
%                 task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
%                passive_stimOn_times, passive_wheel_starts)
% 
%             %% get ecah cell location
%             m2Loc = contains(probe_areas, 'Secondary motor area');
%             m2cells = template_depths > min(min(probe_area_boundaries{m2Loc})) & template_depths < max(max(probe_area_boundaries{m2Loc}));
%             accLoc = contains(probe_areas, 'Anterior cingulate');
%             acccells = template_depths > min(min(probe_area_boundaries{accLoc})) & template_depths < max(max(probe_area_boundaries{accLoc}));
% 
%             %% for each cell, get if significant + psth, plot (1) both sorted by stim response (2) both sorted by mov response (3) scatter of increase post stim and post-move.
%             uniK = unique(spike_templates);
%             Pvalstim = [];
% 
%             for iUnit = 1:length(uniK)
%                 %thisC = ff(im2cells);
%                 thisUnit = uniK(iUnit);
%                 [raster_x, raster_y, t, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, passive_spike_timeline, ...
%                     passive_stimOn_times, ones(length(passive_stimOn_times), 1));
%                 frB = sum(curr_raster_sorted(1:2:end, 100:290), 2);
%                 frA = sum(curr_raster_sorted(1:2:end, 291:481), 2);
%                 p = signrank(frB, frA);
%                 frB = sum(curr_raster_sorted(2:2:end, 100:290), 2);
%                 frA = sum(curr_raster_sorted(2:2:end, 291:481), 2);
%                 p2 = signrank(frB, frA);
%                 Pvalstim = [Pvalstim, p<0.01 & p2<0.01];
%             end
%             uniK_stim = uniK(Pvalstim==1);
%             psthGUI(spike_templates(ismember(spike_templates, uniK_stim)), task_spike_timeline(ismember(spike_templates, uniK_stim)), task_stimOn_times, task_wheel_starts,...
%                 task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline(ismember(spike_templates, uniK_stim)), ...
%                passive_stimOn_times, passive_wheel_starts)
%            
%            
%            
% 
% 
%             uniK = unique(spike_templates);
% 
%             ff = find(acccells);
%             for iacccells = 1:length(find(acccells))
%                 thisC = ff(iacccells);
%                 thisUnit = uniK(iacccells);
%                 [raster_x, raster_y, t, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, passive_spike_timeline, ...
%                     passive_stimOn_times, ones(length(passive_stimOn_times), 1));
%                 accCellPSTHstim = [accCellPSTHstim; curr_smoothed_psth];
%                 frB = sum(curr_raster_sorted(1:1:end, 100:290), 2);
%                 frA = sum(curr_raster_sorted(1:1:end, 391:581), 2);
%                 p = signrank(frB, frA);
%                 accCellPvalstim = [accCellPvalstim, p];
%                 [raster_x, raster_y, t, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, task_spike_timeline, ...
%                     task_wheel_starts, ones(length(task_wheel_starts), 1));
%                 accCellPSTHmove = [accCellPSTHmove; curr_smoothed_psth];
%                 frB = sum(curr_raster_sorted(1:1:end, 100:290), 2);
%                 frA = sum(curr_raster_sorted(1:1:end, 391:581), 2);
%                 p = signrank(frB, frA);
%                 accCellPvalmove = [accCellPvalmove, p];
%             end
% 
%             %% save values for this script, and cell location
%             %get locations (cf script bit andy set on slack)
%             %plot two raster / cell, all waveforms on side in black, this cell in blue
%             %(see if it's more wide or more narrow hehe), ACG same (don't do histogram,
%             %do lines)
%         
%     end
% 
%     %% compare labeling to AP noise manual curation
%     % % QQ check correct matching
%     % AP_load_experimentJF;
%     % APnoiseUnits = good_templates == 0;
%     % removeThese = ~ismember(1:max(spikeTemplates), unique(spikeTemplates));
%     % APnoiseUnits(removeThese) = [];
%     % BCbadUnits = goodUnits == 0;
%     % fracConcordance = sum(BCbadUnits(APnoiseUnits) == 1) / numel(BCbadUnits(APnoiseUnits));
%     %
%     % diffLabeled = find(BCbadUnits(APnoiseUnits) == 1 );
%     % %plot the waveform of units not indentified as noise by me
%     % for iDiffLabeledUnit = 1:length(diffLabeled)
%     %     figure();
%     %     minWv = max([-2, -qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)) + 1]);
%     %     maxWv = min([6-abs(minWv), size(templateWaveforms,3) - qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))]);
%     %     waveformSelect = abs(maxWv)-6:1:maxWv;
%     %     yLim = [min(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)))), ...
%     %         max(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))))];
%     %
%     %     for iSubPlot = 1:6
%     %         subplot(3,2,iSubPlot)
%     %         plot(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, ...
%     %             qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))+waveformSelect(iSubPlot)))
%     %         xlim([0 82])
%     %         ylim([yLim(1), yLim(2)])
%     %         box off; makepretty;
%     %         set(gca,'xtick',[])
%     %         set(gca,'ytick',[])
%     %     end
 end
% 
%     %% look through good/bad cells
% 
%     %% wide vs narrow
end
% 
% smooth_filt = [10,100]; % (trials x frames)
%         
% figure('color', 'white');
% suptitle('m2')
% subplot(221)
% im = nanconv((m2CellPSTHstim(m2CellPvalstim<0.05, :) - nanmean(m2CellPSTHstim(m2CellPvalstim<0.05, 1:200), 2)) ./ nanmean(m2CellPSTHstim(m2CellPvalstim<0.05, 1:200), 2), smooth_filt, 'same');
% [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from stim')
% ylabel('cell # (sorted by stim resp)')
% caxis([-500, 500])
% 
% subplot(222)
% im = nanconv((m2CellPSTHmove(m2CellPvalstim<0.05,:) - nanmean(m2CellPSTHmove(m2CellPvalstim<0.05, 1:200), 2)) ./ nanmean(m2CellPSTHmove(m2CellPvalstim<0.05, 1:200), 2), smooth_filt, 'same');
% %[~, sorted_im_idx] = sort(nanmean(im(:,300:500),2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from move')
% ylabel('cell # (sorted by stim resp)')
% caxis([-500, 500])
% 
% subplot(223)
% im = nanconv((m2CellPSTHmove(m2CellPvalmove<0.05,:) - nanmean(m2CellPSTHmove(m2CellPvalmove<0.05, 1:200), 2)) ./ nanmean(m2CellPSTHmove(m2CellPvalmove<0.05, 1:200), 2), smooth_filt, 'same');
% [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from move')
% ylabel('cell # (sorted by move resp)')
% caxis([-500, 500])
% 
% subplot(224)
% im = nanconv((m2CellPSTHstim(m2CellPvalmove<0.05,:) - nanmean(m2CellPSTHstim(m2CellPvalmove<0.05, 1:200), 2)) ./ nanmean(m2CellPSTHstim(m2CellPvalmove<0.05, 1:200), 2), smooth_filt, 'same');
% %[~, sorted_im_idx] = sort(nanmean(im(:,300:500),2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from stim')
% ylabel('cell # (sorted by move resp)')
% caxis([-500, 500])
% 
% figure('color', 'white');
% suptitle('acc')
% subplot(221)
% im = nanconv((accCellPSTHstim(accCellPvalstim<0.05,:) - nanmean(accCellPSTHstim(accCellPvalstim<0.05, 1:200), 2)) ./ nanmean(accCellPSTHstim(accCellPvalstim<0.05, 1:200), 2), smooth_filt, 'same');
% [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from stim')
% ylabel('cell # (sorted by stim resp)')
% caxis([-500, 500])
% 
% subplot(222)
% im = nanconv((accCellPSTHmove(accCellPvalstim<0.05,:) - nanmean(accCellPSTHmove(accCellPvalstim<0.05, 1:200), 2)) ./ nanmean(accCellPSTHmove(accCellPvalstim<0.05, 1:200), 2), smooth_filt, 'same');
% %[~, sorted_im_idx] = sort(nanmean(im(:,300:500),2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from move')
% ylabel('cell # (sorted by stim resp)')
% caxis([-500, 500])
% 
% subplot(223)
% im = nanconv((accCellPSTHmove(accCellPvalmove<0.05,:) - nanmean(accCellPSTHmove(accCellPvalmove<0.05, 1:200), 2)) ./ nanmean(accCellPSTHmove(accCellPvalmove<0.05, 1:200), 2), smooth_filt, 'same');
% [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from move')
% ylabel('cell # (sorted by move resp)')
% caxis([-500, 500])
% 
% subplot(224)
% im = nanconv((accCellPSTHstim(accCellPvalmove<0.05,:) - nanmean(accCellPSTHstim(accCellPvalmove<0.05, 1:200), 2)) ./ nanmean(accCellPSTHstim(accCellPvalmove<0.05, 1:200), 2), smooth_filt, 'same');
% %[~, sorted_im_idx] = sort(nanmean(im(:,300:500),2));
% colormap(brewermap([], '*RdBu'));
% imagesc(t, [], im(sorted_im_idx, :))
% xlabel('time from stim')
% ylabel('cell # (sorted by move resp)')
% caxis([-500, 500])
% 
% 
% accCellSinc = nanmean((accCellPSTHstim(:, 400:700) - nanmean(accCellPSTHstim(:, 1:200), 2))./nanmean(accCellPSTHstim(:, 1:200), 2), 2);
% accCellMinc = nanmean((accCellPSTHmove(:, 400:700) - nanmean(accCellPSTHmove(:, 1:200), 2))./nanmean(accCellPSTHmove(:, 1:200), 2), 2);
% m2CellSinc = nanmean((m2CellPSTHstim(:, 400:700) - nanmean(m2CellPSTHstim(:, 1:200), 2))./nanmean(m2CellPSTHstim(:, 1:200), 2), 2);
% m2CellMinc = nanmean((m2CellPSTHmove(:, 400:700) - nanmean(m2CellPSTHmove(:, 1:200), 2))./nanmean(m2CellPSTHmove(:, 1:200), 2), 2);
% 
% figure();
% subplot(121)
% title('m2'); hold on;
% scatter(m2CellSinc, m2CellMinc, 2,'filled')
% xlabel('increase after stim')
% ylabel('increase after move')
% xlim([-2, 10])
% ylim([-2, 10])
% makepretty;
% subplot(122)
% title('acc'); hold on; 
% scatter(accCellSinc, accCellMinc, 2,'filled')
% xlabel('increase after stim')
% ylabel('increase after move')
% xlim([-2, 10])
% ylim([-2, 10])
% makepretty;
