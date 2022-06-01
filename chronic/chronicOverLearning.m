
%% download ggl sheet
r=GetGoogleSpreadsheet('1U7mHq17fvQGeBUtLUz2ShpVfvNo63eyk7XRMop5ueJU');

%% animalsAll = {'JF051', 'JF067', 'JF078'};
implantDays = {'2021-08-27', '2022-01-24', '2021-05-20'};
sites = [1, 1, 3];
implantSide = [90, 90, 90];
plotDriftMap = 0;
plotActivity = 1;
ex = [];%[4, 6, 7];
startTrainingDayNum = [NaN, 3, 3];
startStage4DayNum = [NaN, 12, NaN];
% depths1 = [2050, 2800];
% depths2 = [2800, 3500];
depths = {{[0,1200]}, {[2000:3900]}, {[3000,3900;}}; %add 1400:2800: MRN
thisAnimal = 3;

if plotActivity
    bhvData = noGoWorld_behavior(animalsAll{thisAnimal});
end

for iAnimal = thisAnimal%1:size(animalsAll, 2)

    figure(4);
    clf;
    animal = animalsAll{1, iAnimal}; %this animal
    implantDay = implantDays{1, iAnimal}; %this animal
    protocol1 = 'JF_choiceworldStimuli_wheel'; %protocol name contains this name
    protocol2 = 'JF_choiceworldStimuli_wheel_left_center';
    protocol3 = 'nogo';
    protocol0 = 'JF_choiceworldStimuli';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments = AP_find_experimentsJF(animal, protocol0, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol3, flexible_name);
    startTraining = experiments3(1).day;
    startTraining = split(startTraining, '-');
    startTraining = datetime(str2num(startTraining{1}), str2num(startTraining{2}), str2num(startTraining{3}));
    keepPassive = zeros(size(experiments, 1), 1);
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];
    clearvars  alldaysNum singleUnitCount multiUnitCount noiseCount
    if plotActivity
    movingFrac = bhvData(iAnimal).movingFracGo1;
    movingFrac(ex-startTrainingDayNum(iAnimal)+1, :) = [];
    end
    alldays = {experiments.day};
    theseDays = 1:size(experiments, 1);
    theseDays(ex) = [];
    % 29, 30, 31, 
    for curr_day = theseDays
        %curr_day = 35
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
                if isempty(dir([ephysKSfile, filesep, 'sync.mat'])) && plotActivity
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
            else
                multiUnitCount(curr_day) = NaN;
                noiseCount(curr_day) = NaN;
                singleUnitCount(curr_day) = NaN;
                disp('getting quality metrics ...')
                ephysPath = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys', site);
                ephysPath = strrep(ephysPath, 'kilosort2', 'kilosort2');
                ephysap_path = AP_cortexlab_filenameJF(animal, thisDay, experiment, 'ephys_ap', site);
                [spikeTimes, spikeTemplates, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] = ...
                    bc_loadEphysData(ephysPath);
                bc_qualityParamValues; 
                bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
                load(fullfile(savePath, 'qMetric.mat'))
                chronicParamValue; 
                bc_getQualityUnitType;
            end
            

            
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
            loadLFP = 0;
            if plotActivity && ~isempty(dir([ephys_filename, '/site1/spike_templates.npy']))
                AP_load_experimentJF;
            end
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
            if plotActivity && ~isempty(dir([ephys_filename, '/site1/spike_templates.npy']))
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
            if curr_day >= startTrainingDayNum(iAnimal)
                subplot(5, size(experiments, 1), curr_day)
                if curr_day >= startStage4DayNum(iAnimal)
                    hold on;
                    scatter([1, 2, 3], bhvData(iAnimal).goLeft(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2]), [], rgb('DarkBlue'), 'filled')
                    p1 = plot([1, 2, 3], bhvData(iAnimal).goLeft(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2]), 'Color', rgb('DarkBlue'));
                    makepretty;
                    scatter([1, 2, 3], bhvData(iAnimal).noGo(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2]), [], rgb('Red'), 'filled')
                    p2 = plot([1, 2, 3], bhvData(iAnimal).noGo((currDayIdx - startTrainingDayNum(iAnimal) + 1), [1, 3, 2])./ ...
                        bhvData(iAnimal).nTrials(currDayIdx-startTrainingDayNum(iAnimal)+1, [1, 3, 2]), 'Color', rgb('Red'));
                    makepretty;
                    if curr_day == startStage4DayNum(iAnimal)
                        legend([p1, p2], {'\color[rgb]{0,0,1}go left', ' \color[rgb]{1,0,0} no go'})
                    end
                    xticks([1, 2, 3])
                    xticklabels({'Go1', 'NoGo', 'Go2'})
                    ylim([0, 1])


                else
                    plot(bhvData(iAnimal).binBorders(1:end-1), movingFrac(currDayIdx-startTrainingDayNum(iAnimal)+1, :), 'b');
                    box off;
                    if curr_day == startTrainingDayNum(iAnimal)
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
            try
                stimType(trial_conditions(no_move_trials, 1) == 4 & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 1; % go1 stim, contra
                go1Contra = 1;
            catch
                go1Contra = 0;
            end
            try
                stimType(trial_conditions(no_move_trials, 1) == 4 & trial_conditions(no_move_trials, 2) == 0) = 2; % go1 stim, center
                go1Center = 1;
            catch
                go1Center = 0;
            end
            try
                stimType(trial_conditions(no_move_trials, 1) == 7 & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 3; % no go stim, contra
                noGoContra = 1;
            catch
                noGoContra = 0;
            end
            try
                stimType(trial_conditions(no_move_trials, 1) == 7 & trial_conditions(no_move_trials, 2) == 0) = 4; % no go stim, center
                noGoCenter = 1;
            catch
                noGoCenter = 0;
            end
            try
                stimType(ismember(trial_conditions(no_move_trials, 1), [1,2,3,4,5,8:22]) & trial_conditions(no_move_trials, 2) == -implantSide(iAnimal)) = 5; % other stims, contra
                go2Contra = 1;
            catch
                go2Contra = 0;
            end
            try
                stimType(ismember(trial_conditions(no_move_trials, 1), [1,2,3,4,5,8:22]) & trial_conditions(no_move_trials, 2) == 0) = 6; % % other stims, center
                go2Center = 1;
            catch
                go2Center = 0;
            end
            %% center stims top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1))
            colorsO = [rgb('Navy'); rgb('DarkRed'); rgb('Green'); rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            uniquePoss = unique(stimType(~isnan(stimType)));
            spikes1 = spike_depths >= depths1(1) & spike_depths <= depths1(2);
            [~, ~, ~, ~, ~, ~, ~, pltTime4, pltY4, ~] = ...
                plotExCellRasterPSTH(stimType, unique(stimType(~isnan(stimType))), ...
                colorsO, stimOn_times(no_move_trials), spike_templates(spikes1), ...
                spike_times_timeline(spikes1), unique(spike_templates(spikes1)), alphaV, szV, 9, 0, 0);
            pltY4 = (pltY4 - nanmean(pltY4(:, pltTime4 < 0), 2)) ./ nanmean(pltY4(:, pltTime4 < 0), 2);
            
            ll = [1,1,2,2,3,3];
            ll= ll(find(logical(mod(uniquePoss,2))));
            pp = pltY4(logical(mod(uniquePoss,2)), :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(ll(align_group), :)), 1:size(pp(:, :), 1));
            
            if sum(uniquePoss == 1) ==1
                savePSTH(1,curr_day,1, :) = pltY4(uniquePoss==1, :);
            end
            if sum(uniquePoss == 2) ==1
                savePSTH(1,curr_day,2, :) = pltY4(uniquePoss==3, :);
            end
            if sum(uniquePoss == 5) ==1
                savePSTH(1,curr_day,3, :) = pltY4(uniquePoss==5, :);
            end


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
            ll = [1,1,2,2,3,3];
            ll= ll(find(logical(~mod(uniquePoss,2))));
            pp = pltY4(logical(~mod(uniquePoss,2)), :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(ll(align_group), :)), 1:size(pp(:, :), 1));

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
            
            if sum(uniquePoss == 2) ==1
                savePSTH(2,curr_day,1, :) = pltY4(uniquePoss==2, :);
            end
            if sum(uniquePoss == 4) ==1
                savePSTH(2,curr_day,2, :) = pltY4(uniquePoss==4, :);
            end
            if sum(uniquePoss == 6) ==1
                savePSTH(2,curr_day,3, :) = pltY4(uniquePoss==6, :);
            end
            

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
            ll = [1,1,2,2,3,3];
            ll= ll(find(logical(mod(uniquePoss,2))));
            pp = pltY4(logical(mod(uniquePoss,2)), :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(ll(align_group), :)), 1:size(pp(:, :), 1));


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
            
            if sum(uniquePoss == 1) ==1
                savePSTH(3,curr_day,1, :) = pltY4(uniquePoss==1, :);
            end
            if sum(uniquePoss == 3) ==1
                savePSTH(3,curr_day,2, :) = pltY4(uniquePoss==3, :);
            end
            if sum(uniquePoss == 5) ==1
                savePSTH(3,curr_day,3, :) = pltY4(uniquePoss==5, :);
            end
            
            %% contra normalized top striatum
            subplot(5, size(experiments, 1), curr_day+size(experiments, 1)+size(experiments, 1)+size(experiments, 1)+size(experiments, 1))
            colorsO = [rgb('DeepSkyBlue'); ...
                rgb('Crimson'); rgb('MediumSeaGreen')];
            alphaV = 1;
            szV = 2;
            spikes1 = spike_depths >= depths1(1) & spike_depths <= depths1(2);
            %pltY4(2:2:6, :) = (pltY4(2:2:6, :) - nanmean(pltY4(1:2:6, pltTime4 < 0), 2)) ./ nanmean(pltY4(1:2:6, pltTime4 < 0), 2);
            ll = [1,1,2,2,3,3];
            ll= ll(find(logical(~mod(uniquePoss,2))));
            pp = pltY4(logical(~mod(uniquePoss,2)), :);
            psth_lines2 = plot(pltTime4, pp);
            hold on;
            arrayfun(@(align_group) set(psth_lines2(align_group), ...
                'XData', pltTime4, 'YData', pp(align_group, :), ...
                'Color', colorsO(ll(align_group), :)), 1:size(pp(:, :), 1));


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
            
            if sum(uniquePoss == 2) ==1
                savePSTH(4,curr_day,1, :) = pltY4(uniquePoss==2, :);
            end
            if sum(uniquePoss == 4) ==1
                savePSTH(4,curr_day,2, :) = pltY4(uniquePoss==4, :);
            end
            if sum(uniquePoss == 6) ==1
                savePSTH(4,curr_day,3, :) = pltY4(uniquePoss==6, :);
            end
            
            end



        clearvars goodUnits qMetric param
        catch
            singleUnitCount(curr_day) = NaN;
            %deadChannels(curr_day) = NaN;
            multiUnitCount(curr_day) = NaN;
            noiseCount(curr_day) = NaN;
        end
    end
    figure(10+iAnimal);
    clf;
    

    plot(alldaysNum(~isnan(multiUnitCount)), movmean(singleUnitCount(~isnan(singleUnitCount)),5), 'Color', 'k')
    hold on;
    plot(alldaysNum(~isnan(multiUnitCount)),  movmean(multiUnitCount(~isnan(multiUnitCount)),5), 'Color', rgb('Orange'))
    plot(alldaysNum(~isnan(multiUnitCount)),  movmean(noiseCount(~isnan(noiseCount)),15), 'Color', rgb('Red'))
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
    ylim([0 340])
    
    figure(20+iAnimal);
    clf;
    plot(alldaysNum(~isnan(singleUnitCount)), movmean(singleUnitCount(~isnan(singleUnitCount)),5), 'Color', 'k')
    hold on;
    plot(alldaysNum(~isnan(multiUnitCount)),  movmean(multiUnitCount(~isnan(multiUnitCount)),5), 'Color', rgb('Orange'))
    plot(alldaysNum(~isnan(noiseCount)),  movmean(noiseCount(~isnan(noiseCount)),15), 'Color', rgb('Red'))
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
    ylim([0 340])
    


end

figure(); 
averageIncContraDMS = nanmean(savePSTH(1, :, 3, 301:451),4);
averageIncCenterDMS = nanmean(savePSTH(2, :, 3, 301:451),4);

averageIncContraDCS = nanmean(savePSTH(3, :, 3, 301:451),4);
averageIncCenterDCS = nanmean(savePSTH(4, :, 3, 301:451),4);
subplot(411)
scatter(1:26,movmean(averageIncContraDMS([1:13, 15, 17:22, 24:29]), 2), 'filled'); 
subplot(412)
scatter(1:26,movmean(averageIncCenterDMS([1:13, 15, 17:22, 24:29]), 2), 'filled');
subplot(413)
scatter(1:26,movmean(averageIncContraDCS([1:13, 15, 17:22, 24:29]), 2), 'filled');
subplot(414)
scatter(1:26,movmean(averageIncCenterDCS([1:13, 15, 17:22, 24:29]), 2), 'filled');



