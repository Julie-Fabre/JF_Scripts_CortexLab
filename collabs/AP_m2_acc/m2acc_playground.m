
%% load data
animals = {'AP100', 'AP101', 'AP104', 'AP105', 'AP106'};
passiveProtocol = 'AP_lcrGratingPassive';
bhvProtocol = 'AP_stimWheelRight';
theseLocations = {'Secondary motor area', 'Anterior cingulate', 'Prelimbic', 'Infralimbic'};
locationCounts = ones(size(theseLocations, 2), 1);

locPSTHstim = cell(4, 1);
rastersX_stim = cell(4, 1);
rastersY_stim = cell(4, 1);
Pvalstim = cell(4, 1);
locPSTHmove = cell(4, 1);
rastersX_move = cell(4, 1);
rastersY_move = cell(4, 1);
Pvalmove = cell(4, 1);
allWaveforms = cell(4, 1);
allPSS = cell(4, 1);
allACG = cell(4, 1);
allFR = cell(4, 1);
thisCellRec = cell(4, 1); 
recCount = 0; 
for iAnimal = 1:size(animals, 2)
    animal = animals{iAnimal};
    experimentsBhv = AP_find_experimentsJF(animal, bhvProtocol, true);
    experimentsBhv = experimentsBhv([experimentsBhv.ephys]);
    experimentsPass = AP_find_experimentsJF(animal, passiveProtocol, true);
    experimentsPass = experimentsPass([experimentsPass.ephys]);
    for iDay = 1:size(experimentsBhv, 1)
        recCount = recCount + 1;
        day = experimentsBhv(iDay).day;
        experimentBhv = experimentsBhv(iDay).experiment;
        experimentBhvDate = experimentsBhv(iDay).day;
        experimentPassInd = find(cellfun(@(c) ischar(c) && strcmp(c, experimentBhvDate), {experimentsPass(:).day}));
        if ~isempty(experimentPassInd) % we have both passive and task data
            experimentPass = experimentsPass(experimentPassInd).experiment;

            ephysPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys');
            [spikeTimes, spikeTemplates, ...
                templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysPath);
            ephysap_path = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_ap');
            ephysDirPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_dir');
            savePath = fullfile(ephysDirPath, 'qMetrics');

            %% compute quality metrics
            ephysDirPath = AP_cortexlab_filenameJF(animal, day, experimentBhv, 'ephys_dir');
            qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
            rerun = 0;
            if isempty(qMetricsExist) || rerun
                bc_qualityParamValues;
                param.nChannels = 384;
                [qMetric, unitTypes] = bc_runAllQualityMetrics(param, spikeTimes, spikeTemplates, ...
                    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
                load(fullfile(savePath, 'qMetric.mat'))
                load(fullfile(savePath, 'param.mat'))
                clearvars unitType;
                bc_getQualityUnitType;
            else
                load(fullfile(savePath, 'qMetric.mat'))
                load(fullfile(savePath, 'param.mat'))
                clearvars unitType;
                bc_getQualityUnitType;
                %bc_plotGlobalQualityMetric;
            end

            %% look at units responses: passive no movement + bhv spont movement (similar xx)
            clearvars probe_area_boundaries probe_areas
            day = experimentsBhv(iDay).day;
            experiment = experimentBhv;
            isSpikeGlx = 0;
            loadClusters = 0; %don't load AP manually sorted stuff, only bombcell 
            lfp_channel = 'all';
            loadLFP = 1;
            ephys_align = 'cortex';

            AP_load_experiment;
           % AP_get_rewardable_ITI_moves;

            

            task_spike_timeline = spike_times_timeline;
            task_stimOn_times = stimOn_times;
            task_stim_to_move = stim_to_move;
            task_stim_to_feedback = stim_to_feedback;
            task_wheel_starts = move_nostim_rewardable_align; %wheel_starts(wheel_move_iti_idx); %only ITI moves
            task_reward = signals_events.responseTimes(n_trials(1):n_trials(end))';

            %  curr_shank = NaN;
            %  AP_cellrasterJF({stimOn_times}, ...
            %  {trial_conditions(:,1)});

            experiment = experimentPass;
            clearvars spike_times_timeline stimOn_times wheel_starts trial_conditions
            AP_load_experiment;

            passive_spike_timeline = spike_times_timeline;
            passive_stimOn_times = stimOn_times;
            passive_stimIDs = stimIDs;
 %          passive_wheel_starts = wheel_starts;
            passive_wheel_move = wheel_move_time;

            clearvars spike_times_timeline stimOn_times wheel_starts trial_conditions

            % psthGUI(spike_templates, task_spike_timeline, task_stimOn_times, task_wheel_starts, ...
            %    task_wheel_types, task_reward, task_trial_conditions, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
            %    passive_stimOn_times, passive_wheel_starts, ...
            %    passive_wheel_types, passive_trial_conditions)

            %% get each cell's location
            locCells = cell(4, 1);
            for iLocation = 1:size(theseLocations, 2)
                locInd = contains(probe_areas, theseLocations{iLocation});
                if ~isempty(locInd)
                    locCells{iLocation} = template_depths > min(min(probe_area_boundaries{locInd})) & template_depths < max(max(probe_area_boundaries{locInd}));
                end
            end

            %% for each cell, get if significant + psth
            uniqueTemplates = unique(spike_templates);
            for iLocation = 1:size(theseLocations, 2)
                if ~isempty(locCells{iLocation})
                    theseCells = find(locCells{iLocation});
                    for iCell = 1:length(theseCells)
                        thisCellRec{iLocation}(locationCounts(iLocation)) = recCount;
                        thisCell = theseCells(iCell);
                        thisUnit = uniqueTemplates(thisCell);
                        allWaveforms{iLocation}(locationCounts(iLocation), :) = waveforms(thisUnit, :);
                        [ccg, ~] = CCGBz([double(task_spike_timeline(spike_templates == thisUnit)); double(task_spike_timeline(spike_templates == thisUnit))], ...
                            [ones(size((task_spike_timeline(spike_templates == thisUnit)), 1), 1); ...
                            ones(size((task_spike_timeline(spike_templates == thisUnit)), 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function
                        %ephysProp.acg(iUnit, :) = ccg(:, 1, 1);

                        %% compute post spike suppression


                        allPSS{iLocation}(locationCounts(iLocation)) = bc_computePSS(ccg(:, 1, 1));
                        allACG{iLocation}(locationCounts(iLocation), :) = ccg(:, 1, 1);
                        allFR{iLocation}(locationCounts(iLocation)) = sum(spike_templates == thisUnit) ./ (abs(min(task_spike_timeline)-max(task_spike_timeline)));
                        [raster_x, raster_y, t, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, passive_spike_timeline, ...
                            passive_stimOn_times(passive_stimIDs == 3 & isnan(wheel_move_time)), ...
                            ones(length(passive_stimOn_times(passive_stimIDs == 3 & isnan(wheel_move_time))), 1), [-0.3, 0.5], 0.001);
                        locPSTHstim{iLocation} = [locPSTHstim{iLocation}; curr_smoothed_psth];
                        rastersX_stim{iLocation}{locationCounts(iLocation)} = raster_x;
                        rastersY_stim{iLocation}{locationCounts(iLocation)} = raster_y;
                        frB = sum(curr_raster_sorted(1:2:end, 100:290), 2);
                        frA = sum(curr_raster_sorted(1:2:end, 350:540), 2);
                        p = signrank(frB, frA);
                        frB = sum(curr_raster_sorted(2:2:end, 100:290), 2);
                        frA = sum(curr_raster_sorted(2:2:end, 350:540), 2);
                        p2 = signrank(frB, frA);

                        frB = sum(curr_raster_sorted(1:1:end, 100:290), 2);
                        frA = sum(curr_raster_sorted(1:1:end, 350:540), 2);
                        p3 = signrank(frB, frA);

                        %                 clf;
                        %                 subplot(121)
                        %                 plot(t,curr_smoothed_psth)
                        %                 subplot(122)
                        %                 scatter(raster_x, raster_y, 4, 'filled')


                        Pvalstim{iLocation}(locationCounts(iLocation), 1) = p;
                        Pvalstim{iLocation}(locationCounts(iLocation), 2) = p2;
                        Pvalstim{iLocation}(locationCounts(iLocation), 3) = p3;

                        [raster_x, raster_y, t_move, curr_smoothed_psth, trial_sort, curr_raster_sorted] = getRaster(spike_templates, thisUnit, task_spike_timeline, ...
                            task_wheel_starts, ones(length(task_wheel_starts), 1), [-0.7, 0.5], 0.001);
                        locPSTHmove{iLocation} = [locPSTHmove{iLocation}; curr_smoothed_psth];
                        rastersX_move{iLocation}{locationCounts(iLocation)} = raster_x;
                        rastersY_move{iLocation}{locationCounts(iLocation)} = raster_y;
                        frB = sum(curr_raster_sorted(1:2:end, 1:500), 2);
                        frA = sum(curr_raster_sorted(1:2:end, 501:1000), 2);
                        p = signrank(frB, frA);
                        frB = sum(curr_raster_sorted(2:2:end, 1:500), 2);
                        frA = sum(curr_raster_sorted(2:2:end, 501:1000), 2);
                        p2 = signrank(frB, frA);

                        frB = sum(curr_raster_sorted(1:1:end, 1:500), 2);
                        frA = sum(curr_raster_sorted(1:1:end, 501:1000), 2);
                        p3 = signrank(frB, frA);
                        Pvalmove{iLocation}(locationCounts(iLocation), 1) = p;
                        Pvalmove{iLocation}(locationCounts(iLocation), 2) = p2;
                        Pvalmove{iLocation}(locationCounts(iLocation), 3) = p3;
                        locationCounts(iLocation) = locationCounts(iLocation) + 1;
                    end
                end
            end


        end

    end

    %% compare labeling to AP noise manual curation
    % % QQ check correct matching
    % AP_load_experimentJF;
    % APnoiseUnits = good_templates == 0;
    % removeThese = ~ismember(1:max(spikeTemplates), unique(spikeTemplates));
    % APnoiseUnits(removeThese) = [];
    % BCbadUnits = goodUnits == 0;
    % fracConcordance = sum(BCbadUnits(APnoiseUnits) == 1) / numel(BCbadUnits(APnoiseUnits));
    %
    % diffLabeled = find(BCbadUnits(APnoiseUnits) == 1 );
    % %plot the waveform of units not indentified as noise by me
    % for iDiffLabeledUnit = 1:length(diffLabeled)
    %     figure();
    %     minWv = max([-2, -qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)) + 1]);
    %     maxWv = min([6-abs(minWv), size(templateWaveforms,3) - qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))]);
    %     waveformSelect = abs(maxWv)-6:1:maxWv;
    %     yLim = [min(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit)))), ...
    %         max(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))))];
    %
    %     for iSubPlot = 1:6
    %         subplot(3,2,iSubPlot)
    %         plot(templateWaveforms(diffLabeled(iDiffLabeledUnit), :, ...
    %             qMetric.maxChannels(diffLabeled(iDiffLabeledUnit))+waveformSelect(iSubPlot)))
    %         xlim([0 82])
    %         ylim([yLim(1), yLim(2)])
    %         box off; makepretty;
    %         set(gca,'xtick',[])
    %         set(gca,'ytick',[])
    %     end
    % end
end

%% Example cells 
iLocation = 2;
ind = Pvalmove{iLocation}(:, 3) <= 0.01; %m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
find(ind)
%ind = m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;

psthGUI_precalc({rastersX_stim{iLocation}{ind}}, {rastersY_stim{iLocation}{ind}}, locPSTHstim{iLocation}(ind, :), {rastersX_move{iLocation}{ind}}, ...
    {rastersY_move{iLocation}{ind}}, locPSTHmove{iLocation}(ind, :), t, t_move)

%% 75th pcentile example cells
iLocation = 2;
Sinc = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./...
       (nanmean(locPSTHstim{iLocation}(:, 100:290), 2) + nanmean(locPSTHstim{iLocation}(:, 350:540),2)),2);
Minc = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./...
       (nanmean(locPSTHmove{iLocation}(:, 1:500), 2) + locPSTHmove{iLocation}(:, 501:1000)), 2);
PS = prctile(Sinc(Pvalstim{iLocation}(:, 3) <= 0.01 & Pvalmove{iLocation}(:, 3) > 0.01),75);
ind = Pvalstim{iLocation}(:, 3) <= 0.01 & Sinc >= 0.5 & Sinc <= 0.75 & Pvalmove{iLocation}(:, 3) > 0.01; %m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
find(ind)
psthGUI_precalc({rastersX_stim{iLocation}{ind}}, {rastersY_stim{iLocation}{ind}}, locPSTHstim{iLocation}(ind, :), {rastersX_move{iLocation}{ind}}, ...
    {rastersY_move{iLocation}{ind}}, locPSTHmove{iLocation}(ind, :), t, t_move)

PS = prctile(Sinc(Pvalstim{iLocation}(:, 3) <= 0.01 & Pvalmove{iLocation}(:, 3) <= 0.01),75);
PM = prctile(Minc(Pvalstim{iLocation}(:, 3) <= 0.01 & Pvalmove{iLocation}(:, 3) <= 0.01),75);
ind = Pvalstim{iLocation}(:, 3) <= 0.01 & Sinc >= 0.5 & Sinc <= 0.75 & PM >= 0.15 & PM <= 0.30 & Pvalmove{iLocation}(:, 3) <= 0.01; %m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
find(ind)
psthGUI_precalc({rastersX_stim{iLocation}{ind}}, {rastersY_stim{iLocation}{ind}}, locPSTHstim{iLocation}(ind, :), {rastersX_move{iLocation}{ind}}, ...
    {rastersY_move{iLocation}{ind}}, locPSTHmove{iLocation}(ind, :), t, t_move)

PS = prctile(Sinc(Pvalstim{iLocation}(:, 3) <= 0.01),75);
ind = Pvalstim{iLocation}(:, 3) <= 0.01 & Sinc >= 0.5 & Sinc <= 0.75 & Pvalmove{iLocation}(:, 3) > 0.01; %m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
find(ind)
psthGUI_precalc({rastersX_stim{iLocation}{ind}}, {rastersY_stim{iLocation}{ind}}, locPSTHstim{iLocation}(ind, :), {rastersX_move{iLocation}{ind}}, ...
    {rastersY_move{iLocation}{ind}}, locPSTHmove{iLocation}(ind, :), t, t_move)

%% Fraction cells vis/move/both (try bar/line plot) 
figure();
stimAndMoveFraction = nan(11, size(theseLocations, 2));
stimFraction = nan(11, size(theseLocations, 2));
moveFraction = nan(11, size(theseLocations, 2));
for iLocation = 1:size(theseLocations, 2)
    for iRec = unique(thisCellRec{iLocation})
    stimAndMoveFraction(iRec, iLocation) = sum(thisCellRec{iLocation}' == iRec & Pvalstim{iLocation}(:, 3) <= 0.01 & Pvalmove{iLocation}(:, 3) <= 0.01)...
        ./length(Pvalstim{iLocation}(thisCellRec{iLocation} == iRec, 3));
    stimFraction(iRec, iLocation) = sum(thisCellRec{iLocation}' == iRec &Pvalstim{iLocation}(:, 3) <= 0.01 & Pvalmove{iLocation}(:, 3) > 0.01)...
        ./length(Pvalstim{iLocation}(thisCellRec{iLocation} == iRec, 3));
    moveFraction(iRec, iLocation) = sum(thisCellRec{iLocation}' == iRec & Pvalmove{iLocation}(:, 3) <= 0.01 & Pvalstim{iLocation}(:, 3) > 0.01)...
        ./length(Pvalstim{iLocation}(thisCellRec{iLocation} == iRec, 3)); 
    end
end 
plot(nanmean(stimAndMoveFraction), 1:4, 'Color', rgb('Blue') ); hold on;
plot(nanmean(stimFraction), 1:4, 'Color', rgb('Gold') );
plot(nanmean(moveFraction), 1:4, 'Color', rgb('Black') ); 
errorbar(nanmean(stimAndMoveFraction), 1:4, [], [], nanstd(stimAndMoveFraction)./sqrt(size(stimAndMoveFraction,1)), ...
    nanstd(stimAndMoveFraction)./sqrt(size(stimAndMoveFraction,1)),'Color', rgb('Blue'))

errorbar(nanmean(stimFraction), 1:4, [], [], nanstd(stimFraction)./sqrt(size(stimFraction,1)), ...
    nanstd(stimFraction)./sqrt(size(stimFraction,1)),'Color', rgb('Gold'))

errorbar(nanmean(moveFraction), 1:4, [], [], nanstd(moveFraction)./sqrt(size(moveFraction,1)), ...
    nanstd(moveFraction)./sqrt(size(moveFraction,1)),'Color', rgb('Black'))
yticks([1, 2, 3, 4])
yticklabels({theseLocations{1}, theseLocations{2}, theseLocations{3}, theseLocations{4}})
xlabel('fraction of cells')
legend({'stim + move cells', 'stim cells', 'move cells'}, 'Box', 'off')
set(gca, 'YDir', 'Reverse');makepretty;
ylim([0.5, 4.5])
%% Histogram movement responses 
figure();
medThisRecStim =nan( size(theseLocations, 2),11);

medThisRecMove =nan( size(theseLocations, 2),11);

for iLocation = 1:size(theseLocations, 2)
    subplot(size(theseLocations, 2), 1, iLocation)
    % get movement responses for stim cells 
    stimCells = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove =  nanmean((locPSTHmove{iLocation}(stimCells, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells, 1:500), 2) + locPSTHmove{iLocation}(stimCells, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove =  nanmean((locPSTHmove{iLocation}(nonStimCells, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells, 501:1000)), 2);
    % plot
    [countsS, binsS ] = hist(stimCellsMove, -0.8:0.08:0.8);
    h=stairs(binsS, countsS./sum(stimCells), 'LineWidth', 2, 'Color', rgb('MediumSeaGreen')); hold on;
    [countsM, binsM ] = hist(nonStimCellsMove,-0.8:0.08:0.8);
    h1=stairs(binsM, countsM./sum(nonStimCells), 'LineWidth', 2, 'Color', rgb('Black')); hold on;
    xlabel('increase after move onset')
    ylabel('fraction of cells')
    title(theseLocations{iLocation})
    
    makepretty;
    for iRec = unique(thisCellRec{iLocation}) % get median
        % get movement responses for stim cells 
    stimCells = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove =  nanmean((locPSTHmove{iLocation}(stimCells & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCells& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove =  nanmean((locPSTHmove{iLocation}(nonStimCells& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % plot
    [counts, bins ] = hist(stimCellsMove, -0.8:0.08:0.8);
    medThisRecStim(iLocation, iRec) = nanmean(bins(counts == median(counts)));
    [counts, bins ] = hist(nonStimCellsMove,-0.8:0.08:0.8);
    medThisRecMove(iLocation, iRec) = nanmean(bins(counts == median(counts)));
    end
    yy = ylim;
    addV = diff(yy)*0.2;
    line([nanmean(binsS(find(cumsum(countsS) >= floor(median(cumsum(countsS))), 1, 'first')))+ 0.04 ,nanmean(binsS(find(cumsum(countsS) >= floor(median(cumsum(countsS))), 1, 'first'))) + 0.04], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('MediumSeaGreen'), 'LineWidth', 2)
    errorbar(nanmean(binsS(find(cumsum(countsS) >= floor(median(cumsum(countsS))), 1, 'first'))) + 0.04, addV/2+yy(2), [],[],nanstd(medThisRecStim(iLocation,:)),nanstd(medThisRecStim(iLocation,:)),...
        'Color', rgb('MediumSeaGreen'), 'LineWidth',2)
    line([nanmean(binsM(find(cumsum(countsM) >= floor(median(cumsum(countsM))), 1, 'first'))) + 0.04 ,nanmean(binsM(find(cumsum(countsM) >= floor(median(cumsum(countsM))), 1, 'first'))) + 0.04], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('Black'), 'LineWidth', 2)
    errorbar(nanmean(binsM(find(cumsum(countsM) >= floor(median(cumsum(countsM))), 1, 'first'))) + 0.04, addV/2+yy(2), [],[],nanstd(medThisRecMove(iLocation,:)),nanstd(medThisRecMove(iLocation,:)),...
        'Color', rgb('Black'), 'LineWidth',2)
    
    legend({'stim cells', 'non stim cells'}, 'Box', 'off')
end

%% Cum  histogrem 
figure();
medThisRecStim =nan( size(theseLocations, 2),11);

medThisRecMove =nan( size(theseLocations, 2),11);

figure();
    subplot(2, 1, 1)
    iLocation = 1;
    % get movement responses for stim cells 
    stimCells = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove =  nanmean((locPSTHmove{iLocation}(stimCells, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells, 1:500), 2) + locPSTHmove{iLocation}(stimCells, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove =  nanmean((locPSTHmove{iLocation}(nonStimCells, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells, 501:1000)), 2);
    % plot
    iLocation = 2;
    % get movement responses for stim cells 
    stimCells2 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove2 =  nanmean((locPSTHmove{iLocation}(stimCells2, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells2, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells2, 1:500), 2) + locPSTHmove{iLocation}(stimCells2, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells2 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove2 =  nanmean((locPSTHmove{iLocation}(nonStimCells2, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells2, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells2, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells2, 501:1000)), 2);
    % plot
    
    [countsS, binsS ] = hist([stimCellsMove; stimCellsMove2], -0.8:0.08:0.8);
    h=stairs(binsS, countsS./(sum(stimCells) + sum(stimCells2)), 'LineWidth', 2, 'Color', rgb('MediumSeaGreen')); hold on;
    [countsM, binsM ] = hist([nonStimCellsMove; nonStimCellsMove2],-0.8:0.08:0.8);
    h1=stairs(binsM, countsM./(sum(nonStimCells) + sum( nonStimCells2)), 'LineWidth', 2, 'Color', rgb('Black')); hold on;
    xlabel('increase after move onset')
    ylabel('fraction of cells')
    title('M2 + ACC')
        makepretty;
    for iRec = unique(thisCellRec{iLocation}) % get median
        % get movement responses for stim cells 
        iLocation = 1
    stimCellsRec = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec =  nanmean((locPSTHmove{iLocation}(stimCellsRec & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
     iLocation = 2
    stimCellsRec2 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec2 =  nanmean((locPSTHmove{iLocation}(stimCellsRec2 & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec2 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec2 =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
   
    
    % plot
   medThisRecStim(iLocation, iRec) = nanmedian([stimCellsMoveRec;stimCellsMoveRec2]);
    
    medThisRecMove(iLocation, iRec) = nanmedian([nonStimCellsMoveRec; nonStimCellsMoveRec2]);
    end
    yy = ylim;
    addV = diff(yy)*0.3;

    line([nanmedian([stimCellsMove; stimCellsMove2]) ,nanmedian([stimCellsMove; stimCellsMove2]) ], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('MediumSeaGreen'), 'LineWidth', 2)
    %[x y w h]
    
rectangle('Position',[nanmedian([stimCellsMove; stimCellsMove2])-(nanstd(medThisRecStim(iLocation,:))./sqrt(11)),yy(1),2*nanstd(medThisRecStim(iLocation,:))./sqrt(11),yy(2)-yy(1)],...
    'FaceColor',[rgb('MediumSeaGreen'), 0.2],'EdgeColor',[rgb('MediumSeaGreen'), 0],...
    'LineWidth',3)
%     errorbar(nanmedian([stimCellsMove; stimCellsMove2]) , addV/2+yy(2), [],[],nanstd(medThisRecStim(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),nanstd(medThisRecStim(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),...
%         'Color', rgb('MediumSeaGreen'), 'LineWidth',2)
    line([nanmedian([nonStimCellsMove; nonStimCellsMove2])  ,nanmedian([nonStimCellsMove; nonStimCellsMove2]) ], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('Black'), 'LineWidth', 2)
%     errorbar(nanmedian([nonStimCellsMove; nonStimCellsMove2]) , addV/3+yy(2), [],[],nanstd(medThisRecMove(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),nanstd(medThisRecMove(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),...
%         'Color', rgb('Black'), 'LineWidth',2)
    rectangle('Position',[nanmedian([nonStimCellsMove; nonStimCellsMove2])-(nanstd(medThisRecMove(iLocation,:))./sqrt(11)),yy(1),2*nanstd(medThisRecMove(iLocation,:))./sqrt(11),yy(2)-yy(1)],...
    'FaceColor',[rgb('Black'), 0.2],'EdgeColor',[rgb('Black'), 0],...
    'LineWidth',3)
ylim([yy(1), yy(2)])

    
    legend({'stim cells', 'non stim cells'}, 'Box', 'off')
    
    % Shuffle labels

        subplot(2, 1, 2)
    iLocation = 3;
    % get movement responses for stim cells 
    stimCells3 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove3 =  nanmean((locPSTHmove{iLocation}(stimCells3, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells3, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells3, 1:500), 2) + locPSTHmove{iLocation}(stimCells3, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells3 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove3 =  nanmean((locPSTHmove{iLocation}(nonStimCells3, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells3, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells3, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells3, 501:1000)), 2);
    % plot
    iLocation = 4;
    % get movement responses for stim cells 
    stimCells4 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMove4 =  nanmean((locPSTHmove{iLocation}(stimCells4, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCells4, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCells4, 1:500), 2) + locPSTHmove{iLocation}(stimCells4, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCells4 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMove4 =  nanmean((locPSTHmove{iLocation}(nonStimCells4, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCells4, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCells4, 1:500), 2) + locPSTHmove{iLocation}(nonStimCells4, 501:1000)), 2);
    % plot
    
    [countsS, binsS ] = hist([stimCellsMove; stimCellsMove2; stimCellsMove3; stimCellsMove4], -0.8:0.08:0.8);
    h=stairs(binsS, countsS./(sum(stimCells) + sum(stimCells2)+ sum(stimCells3)+ sum(stimCells4)), 'LineWidth', 2, 'Color', rgb('MediumSeaGreen')); hold on;
    [countsM, binsM ] = hist([nonStimCellsMove; nonStimCellsMove2; nonStimCellsMove3; nonStimCellsMove4],-0.8:0.08:0.8);
    h1=stairs(binsM, countsM./(sum(nonStimCells) + sum( nonStimCells2)+ sum( nonStimCells3)+ sum( nonStimCells4)), 'LineWidth', 2, 'Color', rgb('Black')); hold on;
    xlabel('increase after move onset')
    ylabel('fraction of cells')
    title('M2 + ACC + PL + IL')
    for iRec = unique(thisCellRec{iLocation}) % get median
        % get movement responses for stim cells 
        % get movement responses for stim cells 
        iLocation = 1
    stimCellsRec = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec =  nanmean((locPSTHmove{iLocation}(stimCellsRec & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
     iLocation = 2
    stimCellsRec2 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec2 =  nanmean((locPSTHmove{iLocation}(stimCellsRec2 & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec2 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec2 =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec2& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
   
    
        iLocation = 3
    stimCellsRec3 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec3 =  nanmean((locPSTHmove{iLocation}(stimCellsRec3 & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec3& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec3& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec3& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec3 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec3 =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec3& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec3& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec3& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec3& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
     iLocation = 4
    stimCellsRec4 = Pvalstim{iLocation}(:, 3) <= 0.01;
    stimCellsMoveRec4 =  nanmean((locPSTHmove{iLocation}(stimCellsRec4 & thisCellRec{iLocation}' == iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(stimCellsRec4& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(stimCellsRec4& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(stimCellsRec4& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
    % movement responses for non-stim cells 
    nonStimCellsRec4 = Pvalstim{iLocation}(:, 3) > 0.01;
    nonStimCellsMoveRec4 =  nanmean((locPSTHmove{iLocation}(nonStimCellsRec4& thisCellRec{iLocation}'==iRec, 501:1000) - nanmean(locPSTHmove{iLocation}(nonStimCellsRec4& thisCellRec{iLocation}'==iRec, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(nonStimCellsRec4& thisCellRec{iLocation}'==iRec, 1:500), 2) + locPSTHmove{iLocation}(nonStimCellsRec4& thisCellRec{iLocation}'==iRec, 501:1000)), 2);
   
    
    % plot
   medThisRecStim(iLocation, iRec) = nanmedian([stimCellsMoveRec;stimCellsMoveRec2;stimCellsMoveRec3;stimCellsMoveRec4]);
    
    medThisRecMove(iLocation, iRec) = nanmedian([nonStimCellsMoveRec; nonStimCellsMoveRec2;nonStimCellsMoveRec3; nonStimCellsMoveRec4]);
    end
    yy = ylim;
    addV = diff(yy)*0.3;

    line([nanmedian([stimCellsMove; stimCellsMove2;stimCellsMove3; stimCellsMove4]) ,nanmedian([stimCellsMove; stimCellsMove2;stimCellsMove3; stimCellsMove4]) ], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('MediumSeaGreen'), 'LineWidth', 2)
    %errorbar(nanmedian([stimCellsMove; stimCellsMove2;stimCellsMove3; stimCellsMove4]) , addV/2+yy(2), [],[],nanstd(medThisRecStim(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),nanstd(medThisRecStim(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),...
    %    'Color', rgb('MediumSeaGreen'), 'LineWidth',2)
    rectangle('Position',[nanmedian([stimCellsMove; stimCellsMove2;stimCellsMove3; stimCellsMove4])-(nanstd(medThisRecStim(iLocation,:))./sqrt(11)),yy(1),2*nanstd(medThisRecStim(iLocation,:))./sqrt(11),yy(2)-yy(1)],...
    'FaceColor',[rgb('MediumSeaGreen'), 0.2],'EdgeColor',[rgb('MediumSeaGreen'), 0],...
    'LineWidth',3)

    line([nanmedian([nonStimCellsMove; nonStimCellsMove2;nonStimCellsMove3; nonStimCellsMove4])  ,nanmedian([nonStimCellsMove; nonStimCellsMove2;nonStimCellsMove3; nonStimCellsMove4]) ], ...
        [yy(1) addV+yy(2)], 'LineStyle', '--', 'Color', rgb('Black'), 'LineWidth', 2)
    rectangle('Position',[nanmedian([nonStimCellsMove; nonStimCellsMove2;nonStimCellsMove3; nonStimCellsMove4])-(nanstd(medThisRecMove(iLocation,:))./sqrt(11)),yy(1),2*nanstd(medThisRecMove(iLocation,:))./sqrt(11),yy(2)-yy(1)],...
    'FaceColor',[rgb('Black'), 0.2],'EdgeColor',[rgb('Black'), 0],...
    'LineWidth',3)
    %errorbar(nanmedian([nonStimCellsMove; nonStimCellsMove2;nonStimCellsMove3; nonStimCellsMove4]) , addV/3+yy(2), [],[],...
    %nanstd(medThisRecMove(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),nanstd(medThisRecMove(iLocation,:))./sqrt(length(medThisRecStim(iLocation,:))),...
     %   'Color', rgb('Black'), 'LineWidth',2)
     ylim([yy(1), yy(2)])
    makepretty;
    legend({'stim cells', 'non stim cells'}, 'Box', 'off')
%% Global averages
for iLocation = 1:size(theseLocations, 2)
    subplot(211);
    hold on;
    av = nanmean(locPSTHstim{iLocation}(:, :));
    plot(t, av)
    xlabel('time from stim')
    ylabel('mean F.R.')
    xlim([t(1), t(end)])
    makepretty;

    subplot(212);
    hold on;
    av = nanmean(locPSTHmove{iLocation}(:, :));
    plot(t_move, av)
    xlabel('time from move')
    ylabel('mean F.R.')
    xlim([t_move(1), t_move(end)])
    makepretty;

end
legend({theseLocations{1}, theseLocations{2}, theseLocations{3}, theseLocations{4}})
figure();
clf;
for iLocation = 1:2
    templates_max_signfix = bsxfun(@times, allWaveforms{iLocation}, ...
        sign(abs(min(allWaveforms{iLocation}, [], 2))-abs(max(allWaveforms{iLocation}, [], 2))));

    [~, waveform_trough] = min(allWaveforms{iLocation}, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(allWaveforms{iLocation}(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(allWaveforms{iLocation}, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    subplot(2, 2, iLocation);
    hold on;
    %av = nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, :) );
    plot(t, nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, :)))
    plot(t, nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, :)))
    %plot(t, nanmean(locPSTHstim{iLocation}(templateDuration_us > 400 & allPSS{iLocation}' < 50, :) ))
    xlabel('time from stim')
    ylabel('mean F.R.')
    xlim([t(1), t(end)])
    legend({'narrow', 'wide'})

    makepretty;

    subplot(2, 2, iLocation+2);
    hold on;
    plot(t_move, nanmean(locPSTHmove{iLocation}(templateDuration_us <= 400, :)))
    plot(t_move, nanmean(locPSTHmove{iLocation}(templateDuration_us > 400, :)))
    %plot(t_move, nanmean(locPSTHmove{iLocation}(templateDuration_us > 400 & allPSS{iLocation}' < 50, :) ))
    xlabel('time from move')
    ylabel('mean F.R.')
    xlim([t_move(1), t_move(end)])
    makepretty;
    legend({'narrow', 'wide'})

end

%% PSTH all cells stim + move (sorted either way)
figure('Color', 'w');
clf;
suptitle('All cells');
hold on;
for iLocation = 1:size(theseLocations, 2)

    subplot(4, 4, (iLocation - 1)*4+1)
    title([theseLocations{iLocation}, ', ', num2str(size(Pvalstim{iLocation}, 1)), ' cells']);
    hold on;
    %ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHstim{iLocation}(:, :) - nanmean(locPSTHstim{iLocation}(:, 1:200), 2)) ./ (nanmean(locPSTHstim{iLocation}(:, 1:200), 2) + nanmean(locPSTHstim{iLocation}(:, 201:end),2));
    [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], im(sorted_im_idx, :))
    xlabel('time from stim')
    ylabel('cell')
    ylim([0, size(sorted_im_idx, 1) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+2)
    hold on;
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHmove{iLocation}(:, :) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2)) ./ (nanmean(locPSTHmove{iLocation}(:, 1:500), 2) + nanmean(locPSTHmove{iLocation}(:, 501:end),2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], im(sorted_im_idx, :))
    xlabel('time from move')
    ylabel('cell')
    ylim([0, size(sorted_im_idx, 1) + 0.5])
    xlim([t_move(1), t_move(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+3)
    hold on;
    %ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHmove{iLocation}(:, :) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2)) ./ (nanmean(locPSTHmove{iLocation}(:, 1:500), 2) + nanmean(locPSTHmove{iLocation}(:, 501:end),2));

    [~, sorted_im_idx] = sort(nanmean(im(:, 500:1000), 2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], im(sorted_im_idx, :))
    xlabel('time from move')
    ylabel('cell')
    ylim([0, size(sorted_im_idx, 1) + 0.5])
    xlim([t_move(1), t_move(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+4)
    hold on;
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHstim{iLocation}(:, :) - nanmean(locPSTHstim{iLocation}(:, 1:200), 2)) ./ (nanmean(locPSTHstim{iLocation}(:, 1:200), 2) + nanmean(locPSTHstim{iLocation}(:, 201:end),2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], im(sorted_im_idx, :))
    xlabel('time from stim')
    ylabel('cell')
    ylim([0, size(sorted_im_idx, 1) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

end

%% PSTH all cells sorted by waveform
figure('Color', 'w');
clf;
suptitle('All cells, sorted by ''type''');
hold on;
for iLocation = 1:size(theseLocations, 2)

    templates_max_signfix = bsxfun(@times, allWaveforms{iLocation}, ...
        sign(abs(min(allWaveforms{iLocation}, [], 2))-abs(max(allWaveforms{iLocation}, [], 2))));

    [~, waveform_trough] = min(allWaveforms{iLocation}, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(allWaveforms{iLocation}(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(allWaveforms{iLocation}, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;


    subplot(4, 4, (iLocation - 1)*4+1);
    hold on;
    title([theseLocations{iLocation}, ', ', num2str(size(Pvalstim{iLocation}, 1)), ' cells']);
    hold on;
    %ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im1 = (locPSTHstim{iLocation}(templateDuration_us <= 400, :) - nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 1:200), 2)) ...
        ./ (nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 1:200), 2) + nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 201:end),2));
    [~, sorted_im_idx1] = sort(nanmean(im1(:, 300:500), 2));
    im2 = (locPSTHstim{iLocation}(templateDuration_us > 400, :) - nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 1:200), 2)) ...
        ./ (nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 1:200), 2) + nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 201:end),2));
    [~, sorted_im_idx2] = sort(nanmean(im2(:, 300:500), 2));

    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], [im1(sorted_im_idx1, :); im2(sorted_im_idx2, :)])
    xlabel('time from stim')
    ylabel('cell')
    ylim([0, size(sorted_im_idx1, 1) + size(sorted_im_idx2, 1) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+2);
    hold on;
    im1 = (locPSTHmove{iLocation}(templateDuration_us <= 400, :) - nanmean(locPSTHmove{iLocation}(templateDuration_us <= 400, 1:500), 2)) ...
        ./ (nanmean(locPSTHmove{iLocation}(templateDuration_us <= 400, 1:500), 2) + nanmean(locPSTHmove{iLocation}(templateDuration_us <= 400, 501:end),2));
    im2 = (locPSTHmove{iLocation}(templateDuration_us > 400, :) - nanmean(locPSTHmove{iLocation}(templateDuration_us > 400, 1:500), 2)) ...
        ./ (nanmean(locPSTHmove{iLocation}(templateDuration_us > 400, 1:500), 2)  + nanmean(locPSTHmove{iLocation}(templateDuration_us > 400, 501:end),2));

    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], [im1(sorted_im_idx1, :); im2(sorted_im_idx2, :)])
    colormap(brewermap([], '*RdBu'));
    xlabel('time from move')
    ylabel('cell')
    ylim([0, size(sorted_im_idx1, 1) + size(sorted_im_idx2, 1) + 0.5])
    xlim([t_move(1), t_move(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+3)
    hold on;
    [~, sorted_im_idx1] = sort(nanmean(im1(:, 500:1000), 2));
    [~, sorted_im_idx2] = sort(nanmean(im2(:, 500:1000), 2));

    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], [im1(sorted_im_idx1, :); im2(sorted_im_idx2, :)])
    xlabel('time from move')
    ylabel('cell')
    ylim([0, size(sorted_im_idx1, 1) + size(sorted_im_idx2, 1) + 0.5])
    xlim([t_move(1), t_move(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+4);
    hold on;
    im1 = (locPSTHstim{iLocation}(templateDuration_us <= 400, :) - nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 1:200), 2)) ...
        ./ (nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 1:200), 2) + nanmean(locPSTHstim{iLocation}(templateDuration_us <= 400, 201:end),2));
    im2 = (locPSTHstim{iLocation}(templateDuration_us > 400, :) - nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 1:200), 2)) ...
        ./ (nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 1:200), 2) + nanmean(locPSTHstim{iLocation}(templateDuration_us > 400, 201:end),2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], [im1(sorted_im_idx1, :); im2(sorted_im_idx2, :)])
    xlabel('time from stim')
    ylabel('cell')
    ylim([0, size(sorted_im_idx1, 1) + size(sorted_im_idx2, 1) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

end

%% correlation stim/move responses
figure();
for iLocation = 1:size(theseLocations, 2)
    subplot(2,2,iLocation)
    Sinc = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./...
        (nanmean(locPSTHstim{iLocation}(:, 100:290), 2) + nanmean(locPSTHstim{iLocation}(:, 350:540),2)),2);
    Minc = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./...
        (nanmean(locPSTHmove{iLocation}(:, 1:500), 2) + locPSTHmove{iLocation}(:, 501:1000)), 2);
    %scatter(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)))
    [R, P] = corrcoef(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)));
    hist3([Sinc, Minc], {-1:0.2:5, -1:0.2:5}, 'CdataMode', 'auto')
    view([0, 90])
    ylim([-1, 5])
    xlim([-1, 5])
    xlabel('increase after stim')
    ylabel('increase after move')
    title([theseLocations{iLocation}, ' p = ', num2str(P(2))])
    colorbar;
end
%% PSTH and scatter stim/rew
figure('Color', 'w');
clf;
for iLocation = 1:size(theseLocations, 2)

    subplot(4, 4, (iLocation - 1)*4+1)
    title([theseLocations{iLocation}, ', ', num2str(size(Pvalstim{iLocation}, 1)), ' cells']);
    hold on;
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHstim{iLocation}(ind, :) - nanmean(locPSTHstim{iLocation}(ind, 1:200), 2)) ./ nanmean(locPSTHstim{iLocation}(ind, 1:200), 2);
    [~, sorted_im_idx] = sort(nanmean(im(:, 300:500), 2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], im(sorted_im_idx, :))
    xlabel('time from stim')
    ylabel('cell (sorted by stim increase)')
    ylim([0, sum(ind) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+2)
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    im = (locPSTHmove{iLocation}(ind, :) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2)) ./ nanmean(locPSTHmove{iLocation}(ind, 1:500), 2);
    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], im(sorted_im_idx, :))
    xlabel('time from move')
    ylabel('cell (sorted by stim increase)')
    ylim([0, sum(ind) + 0.5])
    xlim([t_move(1), t_move(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+3)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01 & Pvalstim{iLocation}(:, 3) <= 0.01;
    Sinc_sm = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_sm = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    Sinc_s = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_s = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Sinc_m = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_m = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    scatter(Sinc_m, Minc_m, 5, rgb('Red'), 'filled');
    hold on;
    scatter(Sinc_s, Minc_s, 5, rgb('Green'), 'filled');
    scatter(Sinc_sm, Minc_sm, 5, rgb('Purple'), 'filled');
    xlabel('increase after stim')
    ylabel('increase after move')
    xlim([-1, 20])
    ylim([-1, 20])
    %set(gca, 'YScale', 'log')
    %set(gca, 'XScale', 'log')
    legend({'Move cells', 'Stim cells', 'Move+stim cells'})

    subplot(4, 4, (iLocation - 1)*4+4)
    Sinc = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./nanmean(locPSTHstim{iLocation}(:, 100:290), 2), 2);
    Minc = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    %scatter(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)))
    [R, P] = corrcoef(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)));
    hist3([Sinc, Minc], {-1:0.5:20, -1:0.5:20}, 'CdataMode', 'auto')
    view([0, 90])
    title(['P = ', num2str(P(2))])
end

figure('Color', 'w');
clf;
for iLocation = 1:size(theseLocations, 2)

    subplot(4, 4, (iLocation - 1)*4+1)
    title([theseLocations{iLocation}, ', ', num2str(size(Pvalmove{iLocation}, 1)), ' cells']);
    hold on;
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    im = (locPSTHmove{iLocation}(ind, :) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2)) ./ nanmean(locPSTHmove{iLocation}(ind, 1:500), 2);
    [~, sorted_im_idx] = sort(nanmean(im(:, 500:1000), 2));
    colormap(brewermap([], '*RdBu'));
    imagesc(t_move, [], im(sorted_im_idx, :))
    xlabel('time from stim')
    ylabel('cell (sorted by stim increase)')
    ylim([0, sum(ind) + 0.5])
    xlim([t_move(1), t(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+2)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    im = (locPSTHstim{iLocation}(ind, :) - nanmean(locPSTHstim{iLocation}(ind, 1:200), 2)) ./ nanmean(locPSTHstim{iLocation}(ind, 1:200), 2);
    colormap(brewermap([], '*RdBu'));
    imagesc(t, [], im(sorted_im_idx, :))
    xlabel('time from move')
    ylabel('cell (sorted by stim increase)')
    ylim([0, sum(ind) + 0.5])
    xlim([t(1), t(end)])
    colorbar;
    caxis([-4, 4])

    subplot(4, 4, (iLocation - 1)*4+3)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01 & Pvalstim{iLocation}(:, 3) <= 0.01;
    Sinc_sm = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_sm = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    Sinc_s = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_s = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Sinc_m = nanmean((locPSTHstim{iLocation}(ind, 350:540) - nanmean(locPSTHstim{iLocation}(ind, 100:290), 2))./nanmean(locPSTHstim{iLocation}(ind, 100:290), 2), 2);
    Minc_m = nanmean((locPSTHmove{iLocation}(ind, 501:1000) - nanmean(locPSTHmove{iLocation}(ind, 1:500), 2))./nanmean(locPSTHmove{iLocation}(ind, 1:500), 2), 2);
    scatter(Sinc_m, Minc_m, 5, rgb('Red'), 'filled');
    hold on;
    scatter(Sinc_s, Minc_s, 5, rgb('Green'), 'filled');
    scatter(Sinc_sm, Minc_sm, 5, rgb('Purple'), 'filled');
    xlabel('increase after stim')
    ylabel('increase after move')
    xlim([-1, 20])
    ylim([-1, 20])
    %set(gca, 'YScale', 'log')
    %set(gca, 'XScale', 'log')
    legend({'Move cells', 'Stim cells', 'Move+stim cells'})

    subplot(4, 4, (iLocation - 1)*4+4)
    Sinc = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./nanmean(locPSTHstim{iLocation}(:, 100:290), 2), 2);
    Minc = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    %scatter(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)))
    [R, P] = corrcoef(Sinc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)), Minc(~isnan(Sinc) & ~isinf(Sinc) & ~isnan(Minc) & ~isinf(Minc)));
    hist3([Sinc, Minc], {-1:0.5:20, -1:0.5:20}, 'CdataMode', 'auto')
    view([0, 90])
    title(['P = ', num2str(P(2))])
end

figure();
clf;
hold on;

%% wide vs narrow cells
figure('Color', 'w');
figure(8);
clf;
figure(7);
clf;
for iLocation = 1:size(theseLocations, 2)
    figure(7)
    templates_max_signfix = bsxfun(@times, allWaveforms{iLocation}, ...
        sign(abs(min(allWaveforms{iLocation}, [], 2))-abs(max(allWaveforms{iLocation}, [], 2))));

    [~, waveform_trough] = min(allWaveforms{iLocation}, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(allWaveforms{iLocation}(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(allWaveforms{iLocation}, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    subplot(4, 5, (iLocation - 1)*5+1)
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    title([theseLocations{iLocation}, ', ', num2str(sum(ind)), ' stim cells']);
    hold on;


    plot(1:82, allWaveforms{iLocation}(ind & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind & templateDuration_us > 400, :), 'Color', 'b')
    ylabel('a.u.')
    xlabel('time')

    subplot(4, 5, (iLocation - 1)*5+2)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Minc_m = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    title([num2str(sum(ind == 1 & Minc_m > 0)), ' + move cells']);
    hold on;
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m > 0 & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m > 0 & templateDuration_us > 400, :), 'Color', 'b')
    ylabel('a.u.')
    xlabel('time')

    subplot(4, 5, (iLocation - 1)*5+3)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Minc_m = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    title([num2str(sum(ind == 1 & Minc_m < 0)), ' - move cells']);
    hold on;
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m < 0 & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m < 0 & templateDuration_us > 400, :), 'Color', 'b')

    xlabel('time')

    figure(8)
    subplot(4, 2, (iLocation - 1)*2+1)
    Sinc_s = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./nanmean(locPSTHstim{iLocation}(:, 100:290), 2), 2);
    hold on;
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    [N, C] = hist3([templateDuration_us, Sinc_s], {0:50:1000, -2:0.2:20}, 'CdataMode', 'auto');
    density_div = nan(size(N, 1), 1);
    density_div(C{1} <= 400, :) = sum(templateDuration_us <= 400);
    density_div(C{1} > 400, :) = sum(templateDuration_us > 400);
    N = N ./ density_div;
    imagesc(C{1}, C{2}, N')
    c = colorbar;
    c.Label.String = 'neuron density';

    colormap(brewermap([], 'Reds'));

    % %      hist3([templateDuration_us, Sinc_s], {0:50:1000, -1:0.5:20}, 'CdataMode', 'auto')
    % %      view([0, 90])
    %
    %      scatter(templateDuration_us(templateDuration_us <= 400 & ind), Sinc_s(templateDuration_us <= 400 & ind), 4, 'r', 'filled')
    %      scatter(templateDuration_us(templateDuration_us > 400 & ind), Sinc_s(templateDuration_us > 400 & ind), 4, 'b', 'filled')
    xlabel('waveform duration (us)')
    ylabel('fraction change after stim')
    ylim([-2, 10])

    subplot(4, 2, (iLocation - 1)*2+2)
    hold on;
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    %      scatter3(templateDuration_us, Sinc_s, Minc_m, 2, 'filled')
    [N, C] = hist3([templateDuration_us, Minc_m], {0:50:1000, -2:0.2:20}, 'CdataMode', 'auto');
    density_div = nan(size(N, 1), 1);
    density_div(C{1} <= 400, :) = sum(templateDuration_us <= 400);
    density_div(C{1} > 400, :) = sum(templateDuration_us > 400);
    N = N ./ density_div;
    imagesc(C{1}, C{2}, N')
    c = colorbar;
    c.Label.String = 'neuron density';

    colormap(brewermap([], 'Reds'));
    %       %caxis([-10, 10])
    % %       view([45, 45])
    % %      scatter(templateDuration_us(templateDuration_us <= 400 & ind), Minc_m(templateDuration_us <= 400 & ind), 4, 'r', 'filled')
    % %      scatter(templateDuration_us(templateDuration_us > 400 & ind), Minc_m(templateDuration_us > 400 & ind), 4, 'b', 'filled')
    xlabel('waveform duration (us)')
    ylabel('fraction change around move')
    %     %ylabel('fraction increase around move')

    ylim([-2, 10])
    figure(7)

end


for iLocation = 1:2 %size(theseLocations, 2)
    figure(9)
    subplot(2, 1, iLocation)
    templates_max_signfix = bsxfun(@times, allWaveforms{iLocation}, ...
        sign(abs(min(allWaveforms{iLocation}, [], 2))-abs(max(allWaveforms{iLocation}, [], 2))));

    [~, waveform_trough] = min(allWaveforms{iLocation}, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(allWaveforms{iLocation}(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(allWaveforms{iLocation}, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    scatter3(allPSS{iLocation}, templateDuration_us, allFR{iLocation}, 4, 'filled')
    ylabel('trough-to-peak us')
    xlabel('ACG measure (similar to Buzaki)')
    zlabel('FR (spikes/s)')
    title(theseLocations{iLocation})
    zlim([0, 20])
    view([70, 50])
    set(gca, 'xdir', 'reverse')
end

figure('Color', 'w');
figure(8);
clf;
figure(7);
clf;
for iLocation = 1:size(theseLocations, 2)
    figure(7)
    templates_max_signfix = bsxfun(@times, allWaveforms{iLocation}, ...
        sign(abs(min(allWaveforms{iLocation}, [], 2))-abs(max(allWaveforms{iLocation}, [], 2))));

    [~, waveform_trough] = min(allWaveforms{iLocation}, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(allWaveforms{iLocation}(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(allWaveforms{iLocation}, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    subplot(4, 5, (iLocation - 1)*5+1)
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    title([theseLocations{iLocation}, ', ', num2str(sum(ind)), ' stim cells']);
    hold on;


    plot(1:82, allWaveforms{iLocation}(ind & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind & templateDuration_us > 400, :), 'Color', 'b')
    ylabel('a.u.')
    xlabel('time')

    subplot(4, 5, (iLocation - 1)*5+2)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Minc_m = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    title([num2str(sum(ind == 1 & Minc_m > 0)), ' + move cells']);
    hold on;
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m > 0 & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m > 0 & templateDuration_us > 400, :), 'Color', 'b')
    ylabel('a.u.')
    xlabel('time')

    subplot(4, 5, (iLocation - 1)*5+3)
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    Minc_m = nanmean((locPSTHmove{iLocation}(:, 501:1000) - nanmean(locPSTHmove{iLocation}(:, 1:500), 2))./nanmean(locPSTHmove{iLocation}(:, 1:500), 2), 2);
    title([num2str(sum(ind == 1 & Minc_m < 0)), ' - move cells']);
    hold on;
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m < 0 & templateDuration_us <= 400, :), 'Color', 'r')
    plot(1:82, allWaveforms{iLocation}(ind == 1 & Minc_m < 0 & templateDuration_us > 400, :), 'Color', 'b')

    xlabel('time')

    figure(8)
    subplot(4, 2, (iLocation - 1)*2+1)
    Sinc_s = nanmean((locPSTHstim{iLocation}(:, 350:540) - nanmean(locPSTHstim{iLocation}(:, 100:290), 2))./nanmean(locPSTHstim{iLocation}(:, 100:290), 2), 2);
    hold on;
    ind = Pvalstim{iLocation}(:, 3) <= 0.01;
    [N, C] = hist3([allPSS{iLocation}', Sinc_s], {0:5:100, -2:0.2:20}, 'CdataMode', 'auto');
    density_div = nan(size(N, 1), 1);
    density_div(C{1} <= 50, :) = sum(allPSS{iLocation} <= 50);
    density_div(C{1} > 50, :) = sum(allPSS{iLocation} > 50);
    N = N ./ density_div;
    imagesc(C{1}, C{2}, N')
    c = colorbar;
    c.Label.String = 'neuron density';

    colormap(brewermap([], 'Reds'));

    % %      hist3([templateDuration_us, Sinc_s], {0:50:1000, -1:0.5:20}, 'CdataMode', 'auto')
    % %      view([0, 90])
    %
    %      scatter(templateDuration_us(templateDuration_us <= 400 & ind), Sinc_s(templateDuration_us <= 400 & ind), 4, 'r', 'filled')
    %      scatter(templateDuration_us(templateDuration_us > 400 & ind), Sinc_s(templateDuration_us > 400 & ind), 4, 'b', 'filled')
    xlabel('pss (ms)')
    ylabel('fraction change after stim')
    ylim([-2, 10])

    subplot(4, 2, (iLocation - 1)*2+2)
    hold on;
    ind = Pvalmove{iLocation}(:, 3) <= 0.01;
    %      scatter3(templateDuration_us, Sinc_s, Minc_m, 2, 'filled')
    [N, C] = hist3([allPSS{iLocation}', Minc_m], {0:5:10, -2:0.2:20}, 'CdataMode', 'auto');
    density_div = nan(size(N, 1), 1);
    density_div(C{1} <= 50, :) = sum(allPSS{iLocation} <= 50);
    density_div(C{1} > 50, :) = sum(allPSS{iLocation} > 50);
    N = N ./ density_div;
    imagesc(C{1}, C{2}, N')
    c = colorbar;
    c.Label.String = 'neuron density';

    colormap(brewermap([], 'Reds'));
    %       %caxis([-10, 10])
    % %       view([45, 45])
    % %      scatter(templateDuration_us(templateDuration_us <= 400 & ind), Minc_m(templateDuration_us <= 400 & ind), 4, 'r', 'filled')
    % %      scatter(templateDuration_us(templateDuration_us > 400 & ind), Minc_m(templateDuration_us > 400 & ind), 4, 'b', 'filled')
    xlabel('pss (ms)')
    ylabel('fraction change around move')
    %     %ylabel('fraction increase around move')

    ylim([-2, 10])
    figure(7)

end

%% flip through cells
ind = m2CellPvalstim(:, 3) <= 0.01; %m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
find(ind)
%ind = m2CellPvalstim(:,1)<=0.1 & m2CellPvalstim(:,2)<=0.1;
iLocation = 4;
psthGUI_precalc({rastersX_stim{iLocation}{ind}}, {rastersY_stim{iLocation}{ind}}, locPSTHstim{iLocation}(ind, :), {rastersX_move{iLocation}{ind}}, ...
    {rastersY_move{iLocation}{ind}}, locPSTHmove{iLocation}(ind, :), t)