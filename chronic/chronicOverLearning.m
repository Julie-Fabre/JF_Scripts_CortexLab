animalsAll = {'JF067'};
implantDays = {'2022-01-24'};
implantSide = 90;
plotDriftMap = 0;
ex = [4, 6, 7];
startTrainingDayNum = 3;
startStage4DayNum = 12;
depths1 = [2050, 2800];
depths2 = [2800, 3500];

bhvData = noGoWorld_behavior(animalsAll);

for iAnimal = 1:size(animalsAll, 2)

    figure(4);
    clf;
    animal = animalsAll{1, iAnimal}; %this animal
    implantDay = implantDays{1, iAnimal}; %this animal
    protocol1 = 'JF_choiceworldStimuli_wheel'; %protocol name contains this name
    protocol2 = 'JF_choiceworldStimuli_wheel_left_center';
    protocol3 = 'stage';
    protocol0 = 'JF_choiceworldStimuli';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments = AP_find_experimentsJF(animal, protocol0, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol1, flexible_name);
    startTraining = experiments3(1).day;
    startTraining = split(startTraining, '-');
    startTraining = datetime(str2num(startTraining{1}), str2num(startTraining{2}), str2num(startTraining{3}));
    keepPassive = zeros(size(experiments, 1), 1);
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];
    movingFrac = bhvData(iAnimal).movingFracGo1;
    movingFrac(ex-startTrainingDayNum+1, :) = [];
    alldays = {experiments.day};
    theseDays = 1:size(experiments, 1);
    theseDays(ex) = [];
    for curr_day = 1:size(theseDays, 2)
        try
        currDayIdx = theseDays(curr_day);
        alldaysNum(curr_day) = hours((datetime(alldays{currDayIdx}, 'InputFormat', 'yyyy-MM-dd') - datetime(implantDay, 'InputFormat', 'yyyy-MM-dd'))/24);
        % if kilosorted, load the spikes, spike_templates, amplitudes,
        % spike_depths
        thisDay = experiments(currDayIdx).day;
        day = experiments(currDayIdx).day;
        experiment = experiments(currDayIdx).experiment(end);
        thisExperiment = experiments(currDayIdx).experiment(end);
        site = 1;
        isSpikeGlx = 1;
        [ephys_filename, ~] = AP_cortexlab_filenameJF(animal, thisDay, thisExperiment, 'ephys');
        if ~isempty(dir([ephys_filename, '/site1/spike_templates.npy'])) %is has been already kilosorted
            verbose = false; % display load progress and some info figures
            load_parts.cam = false;
            load_parts.imaging = false;
            load_parts.ephys = true;

            site = 1; %1,1; 2,4; 3,7
            if ~isempty(dir([ephys_filename, '/site2/spike_templates.npy']))
                site = 2;
            end
            recording = [];
            loadClusters = 0;
            [ephysAPfile, aa] = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys_ap', site, recording);
            if size(ephysAPfile, 2) == 2 %keep only ap
                ephysAPfile = ephysAPfile{1};
            end
            isSpikeGlx = contains(ephysAPfile, '_g');
            if isSpikeGlx
                [ephysKSfile, ~] = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys', site, recording);
                if isempty(dir([ephysKSfile, filesep, 'sync.mat']))
                    sync = syncFT(ephysAPfile, 385, ephysKSfile);
                end
            end
            ephysDirPath = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys_dir', site);
            savePath = fullfile(ephysDirPath, 'qMetrics');
            qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
            if ~isempty(qMetricsExist)

                load(fullfile(savePath, 'qMetric.mat'))
                ephysap_path = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys_ap', site);
                chronicParamValue; 
                bc_getQualityUnitType;
            end
            
%             uniK = unique(spike_templates);
%             for iUnit = 1:length(uniK)
%                 thisUnit = uniK(iUnit);
%                 [raster_x, raster_y, t, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, spike_times_timeline, ...
%                     stimOn_times, ones(length(passive_stimOn_times), 1));
%                 frB = sum(curr_raster_sorted(1:1:end, 100:290), 2);
%                 frA = sum(curr_raster_sorted(1:1:end, 391:581), 2);
%                 p(iUnit) = signrank(frB, frA);
%             end
%             keepTheseUnits = 
            
            multiUnitCount(curr_day) = length(find(unitType == 2));
            noiseCount(curr_day) = length(find(unitType == 0));
            singleUnitCount(curr_day) = length(find(unitType==1));
            %deadChannels(curr_day) = length(unique(channel_map));
            clearvars unitType %for now, keep multi unit for average 
            unitType = nan(length(qMetric.percSpikesMissing), 1);
            unitType(qMetric.nPeaks > param.maxNPeaks | qMetric.nTroughs > param.maxNTroughs | qMetric.somatic ~= param.somatic ...
                | qMetric.spatialDecaySlope <=  param.minSpatialDecaySlope | qMetric.waveformDuration < param.minWvDuration |...
                qMetric.waveformDuration > param.maxWvDuration  | qMetric.waveformBaseline >= param.maxWvBaselineFraction) = 0; %NOISE or NON-SOMATIC
            unitType(isnan(unitType)') = 1; %SINGLE SEXY UNIT
    
            AP_load_experimentJF;
            if plotDriftMap

                %% drift map
                figure(1);
                subplot(1, size(experiments, 1), curr_day)
                [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap([ephys_filename, '/site1/']);
                plotDriftmap(spikeTimes, template_amplitudes, spikeDepths);
                makeprettyLarge;
                xlim([0, max(spikeTimes)])
                ylim([0, 2880])

                xticks([0, max(spikeTimes)])
                xticklabels({'0', num2str(max(spikeTimes)/60)})
                if curr_day == 1
                    ylabel('Depth from tip (\mum)');
                    xlabel('time (min)')
                    makeprettyLarge;
                else
                    set(gca, 'ytick', [])
                    set(gca, 'yticklabel', [])
                    ylabel('')
                end
            end

            %% units: depths * normalized log rate
            figure(3)
            subplot(1, size(experiments, 1), curr_day)

            norm_spike_n = mat2gray(log10(accumarray(spike_templates, 1)+1));
            unit_dots = plot(norm_spike_n, template_depths, '.k', 'MarkerSize', 20);
            xlim([-0.1, 1]);
            ylim([-10, 3880 + 50]);
            if curr_day == 1
                ylabel('Depth (\mum)')
                xlabel('Normalized log rate')
            else
                set(gca, 'ytick', [])
                set(gca, 'yticklabel', [])
            end
            set(gca, 'Ydir', 'reverse')
            makeprettyLarge;

            %% behavior
            figure(4);
            if curr_day >= startTrainingDayNum
                subplot(5, size(experiments, 1), curr_day)
                if curr_day >= startStage4DayNum
                    hold on;
                    scatter([1, 2, 3], bhvData(iAnimal).goLeft(currDayIdx-startTrainingDayNum+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum+1, [1, 3, 2]), [], rgb('DarkBlue'), 'filled')
                    p1 = plot([1, 2, 3], bhvData(iAnimal).goLeft(currDayIdx-startTrainingDayNum+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum+1, [1, 3, 2]), 'Color', rgb('DarkBlue'));
                    makepretty;
                    scatter([1, 2, 3], bhvData(iAnimal).noGo(currDayIdx-startTrainingDayNum+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum+1, [1, 3, 2]), [], rgb('Red'), 'filled')
                    p2 = plot([1, 2, 3], bhvData(iAnimal).noGo((currDayIdx - startTrainingDayNum + 1), [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum+1, [1, 3, 2]), 'Color', rgb('Red'));
                    makepretty;
                    if curr_day == startStage4DayNum
                        legend([p1, p2], {'\color[rgb]{0,0,1}go left', ' \color[rgb]{1,0,0} no go'})
                    end
                    xticks([1, 2, 3])
                    xticklabels({'Go1', 'NoGo', 'Go2'})
                    ylim([0, 1])


                else
                    plot(bhvData(iAnimal).binBorders(1:end-1), movingFrac(currDayIdx-startTrainingDayNum+1, :), 'b');
                    box off;
                    if curr_day == startTrainingDayNum
                        ylabel('frac. mov.')
                        %legend({'go1, contra','go1, center','other, contra','other, center'})
                    else
                        set(gca, 'ytick', [], 'yticklabel', [], 'XColor', 'white', 'YColor', 'white')
                    end
                    line([0, 0], [0, 1], 'Color', 'k')
                    makepretty;
                    ylim([0, 1])
                    xlim([-5, 10])
                end

            end
            stimType = nan(size(trial_conditions(no_move_trials), 1), 1);
            stimType(trial_conditions(no_move_trials, 1) == 4 & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 1; % go1 stim, contra
            stimType(trial_conditions(no_move_trials, 1) == 4 & trial_conditions(no_move_trials, 2) == 0) = 2; % go1 stim, center
            stimType(trial_conditions(no_move_trials, 1) == 7 & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 3; % no go stim, contra
            stimType(trial_conditions(no_move_trials, 1) == 7 & trial_conditions(no_move_trials, 2) == 0) = 4; % no go stim, center
            stimType(ismember(trial_conditions(no_move_trials, 1), [1, 2, 3, 5, 6]) & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 5; % other stims, contra
            stimType(ismember(trial_conditions(no_move_trials, 1), [1, 2, 3, 5, 6]) & trial_conditions(no_move_trials, 2) == 0) = 6; % % other stims, center

            %% center stims top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1))
            colorsO = [rgb('Navy'); rgb('DarkRed'); rgb('Green'); rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            spikes1 = spike_depths >= depths1(1) & spike_depths <= depths1(2);
            [~, ~, ~, ~, ~, ~, ~, pltTime4, pltY4, ~] = ...
                plotExCellRasterPSTH(stimType, unique(stimType(~isnan(stimType))), ...
                colorsO, stimOn_times(no_move_trials), spike_templates(spikes1), ...
                spike_times_timeline(spikes1), unique(spike_templates(spikes1)), alphaV, szV, 9, 0, 0);
            pltY4 = (pltY4 - nanmean(pltY4(:, pltTime4 < 0), 2)) ./ nanmean(pltY4(:, pltTime4 < 0), 2);
            pp = pltY4(2:2:6, :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(align_group, :)), 1:size(pp(1:2, :), 1));


            xlim([pltTime4(1), 0.25])
            line([0, 0], [-0.37, 1.6], 'Color', 'k')
            box off;

            if curr_day == 1
                legend({'go1, center', 'no go, center', 'other, center'})
            else
                set(gca, 'ytick', [], 'yticklabel', [], 'XColor', 'white', 'YColor', 'white')
            end
            if curr_day == min(2, size(experiments, 1))
                xlabel('time from stim onset')
            end
            makepretty;
            ylim([-0.37, 8])

            %% contra normalized top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1)+size(experiments, 1))
            colorsO = [rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            %pltY4(2:2:6, :) = (pltY4(2:2:6, :) - nanmean(pltY4(1:2:6, pltTime4 > 0), 2)) ./ nanmean(pltY4(1:2:6, pltTime4 > 0), 2);
            pp = pltY4(1:2:6, :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(align_group, :)), 1:size(pp(1:2, :), 1));
            xlim([pltTime4(1), 0.25])
            line([0, 0], [-0.37, 1.6], 'Color', 'k')
            box off;

            if curr_day == 1
                legend({'go1, contra', 'no go, contra', 'other, contra'})
            else
                set(gca, 'ytick', [], 'yticklabel', [], 'XColor', 'white', 'YColor', 'white')
            end
            if curr_day == min(2, size(experiments, 1))
                xlabel('time from stim onset')
            end
            makepretty;
            ylim([-0.37, 4])

            %% center stims top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1)+size(experiments, 1)+size(experiments, 1))
            colorsO = [rgb('Navy'); rgb('DarkRed'); rgb('Green'); rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            spikes2 = spike_depths >= depths2(1) & spike_depths <= depths2(2);
            [~, ~, ~, ~, ~, ~, ~, pltTime4, pltY4, ~] = ...
                plotExCellRasterPSTH(stimType, unique(stimType(~isnan(stimType))), ...
                colorsO, stimOn_times(no_move_trials), spike_templates(spikes2), ...
                spike_times_timeline(spikes2), unique(spike_templates(spikes2)), alphaV, szV, 9, 0, 0);
            pltY4 = (pltY4 - nanmean(pltY4(:, pltTime4 < 0), 2)) ./ nanmean(pltY4(:, pltTime4 < 0), 2);
            pp = pltY4(2:2:6, :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(align_group, :)), 1:size(pp(1:2, :), 1));


            xlim([pltTime4(1), 0.25])
            line([0, 0], [-0.37, 1.6], 'Color', 'k')
            box off;

            if curr_day == 1
                legend({'go1, center', 'no go, center', 'other, center'})
            else
                set(gca, 'ytick', [], 'yticklabel', [], 'XColor', 'white', 'YColor', 'white')
            end
            if curr_day == min(2, size(experiments, 1))
                xlabel('time from stim onset')
            end
            makepretty;
            ylim([-0.37, 4])

            %% contra normalized top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1)+size(experiments, 1)+size(experiments, 1)+size(experiments, 1))
            colorsO = [rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            spikes1 = spike_depths >= depths1(1) & spike_depths <= depths1(2);
            %pltY4(2:2:6, :) = (pltY4(2:2:6, :) - nanmean(pltY4(1:2:6, pltTime4 < 0), 2)) ./ nanmean(pltY4(1:2:6, pltTime4 < 0), 2);
            pp = pltY4(1:2:6, :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(align_group, :)), 1:size(pp(1:2, :), 1));
            xlim([pltTime4(1), 0.25])
            line([0, 0], [-0.37, 1.6], 'Color', 'k')
            box off;

            if curr_day == 1
                legend({'go1, contra', 'no go, contra', 'other, contra'})
            else
                set(gca, 'ytick', [], 'yticklabel', [], 'XColor', 'white', 'YColor', 'white')
            end
            if curr_day == min(2, size(experiments, 1))
                xlabel('time from stim onset')
            end
            makepretty;
            ylim([-0.37, 4])


        else
            singleUnitCount(curr_day) = NaN;
            %deadChannels(curr_day) = NaN;
            multiUnitCount(curr_day) = NaN;
            noiseCount(curr_day) = NaN;
        end
        clearvars goodUnits qMetric param
        catch
            singleUnitCount(curr_day) = NaN;
            %deadChannels(curr_day) = NaN;
            multiUnitCount(curr_day) = NaN;
            noiseCount(curr_day) = NaN;
        end
    end
    figure(2);
    clf;
    plot(alldaysNum, singleUnitCount, 'Color', 'k')
    hold on;
    plot(alldaysNum, multiUnitCount, 'Color', rgb('Orange'))
    plot(alldaysNum, noiseCount, 'Color', rgb('Red'))
    grid minor;
    ylabel('# of units')
    makepretty;
    legend({'nice single units', 'multi units', 'noise units'})
    %yyaxis right;
    %plot(alldaysNum, 384-deadChannels, 'Color', 'b')
    %ax = gca;
    %ax.YAxis(2).Color = 'b';
    xlabel('post-implant day #')
    %ylabel('# of channels with no units')
    xlim([alldaysNum(1), alldaysNum(end)])


end