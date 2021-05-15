clear all;
clc;
% to add first half, second half
% add weight %
% 2021 04 08 added no go analysis
animalsAll = {'JF036','JF037','JF038','JF039','JF042', 'JF043', 'JF044'}; 
for iAnimal = 5:size(animalsAll, 2)
    
    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'location'; %protocol name contains this name
    protocol2 = 'goNoGo';
    protocol3 = 'Hack';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments1 = AP_find_experimentsJF(animal, protocol, flexible_name);
    experiments2 = AP_find_experimentsJF(animal, protocol2, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol3, flexible_name);
    experiments(1:size(experiments1, 1)) = experiments1;
    experiments(size(experiments1, 1)+1:size(experiments1, 1)+size(experiments2, 1)) = experiments2;
    experiments(size(experiments1, 1)+size(experiments2, 1)+1:size(experiments1, 1)+size(experiments2, 1)+size(experiments3, 1)) = experiments3;
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];
    if iAnimal == 4
        theseD = [1:19, 21:length(experiments)];
    else
        theseD = 1:length(experiments);
    end
    for curr_day = theseD%[1:19, 21:length(experiments)]%JF039
        
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        trNum = [];
        % If multiple experiments, only use the largest one (usually multiple
        % happens if mess ups first/last one is good)
        %if curr_day ~= length(experiments)
        for curr_experiment = 1:length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
            if ~isempty(block_filename)
                load(block_filename)
                
                trNum(curr_experiment) = length(block.events.newTrialTimes);
            else
                trNum(curr_experiment) = 0;
            end
            
        end
        %
        for curr_experiment = find(trNum == max(trNum))
            
            experiment = experiment_num(curr_experiment);
            
            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
            load(block_filename)
            %correct non repeat trials
            response_trials = 1:length(block.events.responseValues);
            if isfield(block.events, 'repeatTrialValues')
                repeatOnMisses = block.events.repeatTrialValues(response_trials) > 0;
            else
                repeatOnMisses = block.events.repeatNumValues(response_trials) > 0;
            end
            
            %correctTrialsLeft = block.events.trialSideValues(response_trials) == -1 & block.events.hitValues(response_trials) == 1;
            if isfield(block.events, 'trialSideValues')
                correctTrialsRight = block.events.trialSideValues(response_trials) == 1 & block.events.hitValues(response_trials) == 1;
            else
                correctTrialsRight = block.events.feedbackValues(response_trials) == 1; %go1 and go2
            end
            sum(correctTrialsRight);
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);
            
            
            response_trials = 1:length(block.events.endTrialValues);
            clearvars correctTrials
            iImg = 1;
            if isfield(block.events, 'trialSideValues')
                %leftGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == -1));
                rightGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == 1));
                %leftT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == -1));
                rightT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == 1));
            elseif isfield(block.events, 'noGoQuiescenceTimes')
                block.events.trialSideValues(response_trials) = 1;
                hitTimesPossib = [block.events.feedbackTimes, block.events.noGoQuiescenceTimes];
                hitValuesUnsorted = [block.events.feedbackValues, 0 * ones(size(block.events.noGoQuiescenceTimes, 2), 1)'];
                [sortedV, sortedIdx] = sort(hitTimesPossib);
                block.events.hitValues(response_trials) = hitValuesUnsorted(sortedIdx(response_trials));
            else
                block.events.trialSideValues(response_trials) = 1;
                ff = block.events.responseValues;
                
                block.events.hitValues(response_trials) = ff(response_trials);
            end
            if contains(block.expDef, 'goNoGo')
                go1Trials = block.events.trialTypeValues(response_trials) == 1;
                go2Trials = block.events.trialTypeValues(response_trials) == 2;
                noGoTrials = block.events.trialTypeValues(response_trials) == 3;
            elseif contains(block.expDef, 'Hack')
                go1Trials = block.events.stimulusTypeValues(response_trials) == 1;
                go2Trials = block.events.stimulusTypeValues(response_trials) == 2;
                noGoTrials = block.events.stimulusTypeValues(response_trials) == 3;
            else
                go1Trials = ones(size(block.events.trialSideValues(response_trials), 2), 1)';
                go2Trials = zeros(size(block.events.trialSideValues(response_trials), 2), 1)';
                noGoTrials = zeros(size(block.events.trialSideValues(response_trials), 2), 1)';
            end
            if isfield(block.events, 'repeatTrialValues')
                repeatOnMisses = block.events.repeatTrialValues(response_trials) > 0;
            else
                repeatOnMisses = block.events.repeatNumValues(response_trials) > 1;
            end
            correctTrials((iImg - 1)*3+1) = sum(block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & go1Trials ...
                ); % go 1
            correctTrials((iImg - 1)*3+2) = sum(block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & go2Trials ...
                ); % go 2
            correctTrials((iImg - 1)*3+3) = sum(block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & noGoTrials ...
                ); % no go
            
            nTrials((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & go1Trials ...
                );
            nTrials((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & go2Trials ...
                );
            nTrials((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & noGoTrials ...
                );
            
            conditi = [1, 2, 3];
            
            
            goLeft((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & go1Trials ...
                );
            goLeft((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & go2Trials ...
                );
            goLeft((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 & noGoTrials ...
                );
            
            
            if contains(block.expDef, 'goNoGo')
                allOffs = sort([block.events.goStimOffTimes, block.events.noGoStimOffTimes]);
                timeThis = allOffs(response_trials) - block.events.goStimOnTimes(response_trials);
            elseif isfield(block.events, 'stimulusOnTimes')
                timeThis = block.events.feedbackTimes(response_trials) - block.events.stimulusOnTimes(response_trials);
            else
                timeThis = block.events.stimOffTimes(response_trials) - block.events.stimOnTimes(response_trials);
            end
            trialTime((iImg - 1)*3+1) = nanmean(timeThis(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                go1Trials ...
                ));
            trialTime((iImg - 1)*3+2) = nanmean(timeThis(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                go2Trials ...
                ));
            trialTime((iImg - 1)*3+3) = nanmean(timeThis(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                noGoTrials ...
                ));
            
            
            
            noGo((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                (block.events.hitValues(response_trials) == 0 ) & go1Trials ...
                );
            noGo((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                (block.events.hitValues(response_trials) == 0 ) & go2Trials ...
                );
            noGo((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                (block.events.hitValues(response_trials) == 0 ) & noGoTrials ...
                );

            
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1 / wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, wheel_t_resample);
            
            wheel_smooth_t = 0.05; % seconds
            wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
            wheel_velocity = interp1(conv(wheel_t_resample, [1, 1]/2, 'valid'), ...
                diff(smooth(wheel_values_resample, wheel_smooth_samples)), wheel_t_resample)';
            wheel_thresh = 0.025;
            
            timeBin = [-find(wheel_t_resample > 3, 1) - find(wheel_t_resample > 1, 1), find(wheel_t_resample > 3, 1) - find(wheel_t_resample > 1, 1)];
            
            wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
                abs(wheel_velocity(2:end)) > wheel_thresh);
            
            movingRN = abs(wheel_velocity) > 0.022;
            %             figure();
            %             plot(movingRN)
            %             hold on;
            %             plot(wheel_velocity)
            
            movingRN_times = unique(wheel_t_resample(movingRN));
            if isfield(block.events, 'stimOnTimes')
                win = [-5, 20];
                binSize = 0.1;
                binBorders = win(1):binSize:win(2);
                binArray = [];
                for r = 1:length(block.events.stimOnTimes) - 1
                    [n, binCenters] = histdiff(movingRN_times, block.events.stimOnTimes(r), binBorders);
                    binArray(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp = cumsum(binArray > 1);
                movingFrac = cp(end, :) / size(binArray, 1);
            elseif isfield(block.events, 'noGoStimOnTimes')
                win = [-5, 20];
                binSize = 0.1;
                binBorders = win(1):binSize:win(2);
                binArrayGo1 = [];
                theseT = block.events.noGoStimOnTimes(go1Trials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayGo1(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp = cumsum(binArrayGo1 > 1);
                movingFracGo1 = cp(end, :) / size(binArrayGo1, 1);
                
                binArrayGo2 = [];
                theseT = block.events.noGoStimOnTimes(go2Trials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayGo2(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp = cumsum(binArrayGo2 > 1);
                movingFracGo2 = cp(end, :) / size(binArrayGo2, 1);
                
                binArrayNoGo = [];
                theseT = block.events.noGoStimOnTimes(noGoTrials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayNoGo(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp2 = cumsum(binArrayNoGo > 1);
                movingFracNoGo = cp2(end, :) / size(binArrayNoGo, 1);
            else
                win = [-5, 20];
                binSize = 0.1;
                binBorders = win(1):binSize:win(2);
                binArrayGo1 = [];
                theseT = block.events.stimulusOnTimes(go1Trials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayGo1(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp = cumsum(binArrayGo1 > 1);
                movingFracGo1 = cp(end, :) / size(binArrayGo1, 1);
                
                if sum(go2Trials)>=1
                binArrayGo2 = [];
                theseT = block.events.stimulusOnTimes(go2Trials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayGo2(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                
                cp = cumsum(binArrayGo2 > 1);
                movingFracGo2 = cp(end, :) / size(binArrayGo2, 1);
                end
                binArrayNoGo = [];
                theseT = block.events.stimulusOnTimes(noGoTrials);
                for r = 1:length(theseT) - 1
                    [n, binCenters] = histdiff(movingRN_times, theseT(r), binBorders);
                    binArrayNoGo(r, :) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
                end
                cp2 = cumsum(binArrayNoGo > 1);
                if length(theseT) > 1
                    movingFracNoGo = cp2(end, :) / size(binArrayNoGo, 1);
                else
                    movingFracNoGo = nan(1,250);
                end
            end
            
            
            bhv.correctTrials(curr_day, :) = correctTrials;
            bhv.nTrials(curr_day, :) = nTrials;
            bhv.n_trials(curr_day, :) = n_trials;
            bhv.total_water(curr_day, :) = total_water;
            bhv.conditions(curr_day, :) = conditi;
            bhv.goLeft(curr_day, :) = goLeft;
            bhv.noGo(curr_day, :) = noGo;
            bhv.trialTime(curr_day, :) = trialTime;
            keep_day = [keep_day, curr_day];
            if contains(block.expDef, 'NoGo')
                noGoDay = [noGoDay, curr_day];
                if sum(go2Trials) >= 1
                bhv.movingFracGo2(curr_day, :) = movingFracGo2;
                end
                bhv.movingFracGo1(curr_day, :) = movingFracGo1;
                bhv.movingFracNoGo(curr_day, :) = movingFracNoGo;
            else
                bhv.movingFrac(curr_day,:) = movingFrac;
            end
            
            %stims
        end
    end
    
    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels = cellfun(@(day) [' ', day(6:end)], ...
        {experiments(keep_day).day}, 'uni', false);
%     
    figure('Name', animal);
    subplot(321)
    yyaxis left
    
    scatter(day_num, bhv.n_trials(keep_day, :));
    hold on;
    plot(day_num, bhv.n_trials(keep_day, :), 'linewidth', 2);
    ylabel('Trials');
    yyaxis right
    
    scatter(day_num, bhv.total_water(keep_day, :), 'linewidth', 2);
    hold on;
    plot(day_num, bhv.total_water(keep_day, :), 'linewidth', 2);
    ax = gca;
    hold on;
    
    ylabel('Total water (ul)');
    xlabel('Session');
%     set(gca, 'XTick', day_num);
%     set(gca, 'XTickLabel', day_labels);
%     set(gca, 'XTickLabelRotation', 90);
    makepretty;
    
    subplot(323)
    con = bhv.conditions(keep_day(1), :);
    ccc = num2str(con(1, :));
    cc = textscan(ccc, '%s', 'Delimiter', ' ');
    cc = cc{1};
    cc(strcmp('', cc)) = [];
    im = imagesc(1:length(con), 1:size(bhv.correctTrials(keep_day, :), 1), bhv.goLeft(keep_day, :)./bhv.nTrials(keep_day, :));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'Go left (frac)');
    xlabel('image type');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);
    %
    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.correctTrials(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    %axis square;
    caxis([0, 1])
    makepretty;
    
    subplot(324)
    con = bhv.conditions(keep_day(1), :);
    ccc = num2str(con(1, :));
    cc = textscan(ccc, '%s', 'Delimiter', ' ');
    cc = cc{1};
    cc(strcmp('', cc)) = [];
    im = imagesc(1:length(con), 1:size(bhv.correctTrials(keep_day, :), 1), bhv.noGo(keep_day, :)./bhv.nTrials(keep_day, :));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'No go (frac)');
    xlabel('image type');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);
    
    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.noGo(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    %axis square;
    caxis([0, 1])
    makepretty;
    
    subplot(322)
    im = imagesc(1:length(con), 1:size(bhv.trialTime(keep_day, :), 1), bhv.trialTime(keep_day, :));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'Trial time (dirty proxy for reaction time)');
    xlabel('image type');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);
    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.trialTime(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    makepretty;
    
    subplot(325)
    
    transparencyValues = 0:1 / length(noGoDay):1;
   
    
    if length(noGoDay)==0
        transparencyValues = 0:1 / length(keep_day):1;
        for iDay = 1:length(keep_day)
        p1 = plot(binBorders(1:end-1), bhv.movingFrac(keep_day(iDay), :), 'b');
        p1.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
        end
        ylim([0 1])
    else
         for iDay = 1:length(noGoDay)
        p1 = plot(binBorders(1:end-1), bhv.movingFracGo1(noGoDay(iDay), :), 'b');
        p1.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
        ylim([0 1])
        
         end
         if isfield(bhv, 'movingFracGo2')
    p2 = plot(binBorders(1:end-1), bhv.movingFracGo2(noGoDay(iDay), :), 'r');
        p2.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
         legend([p1, p2], {'Go1', 'Go2'})
         end
         ylim([0 1])
    end
    xlim([binBorders(1), binBorders(end)])
   
    xlabel('time from stim onset (s)')
    ylabel('fraction moving')
    makepretty;
    
    subplot(326)
    if length(noGoDay)==0
        else
    transparencyValues = 0:1 / length(noGoDay):1;
    for iDay = 1:length(noGoDay)
        
        p3 = plot(binBorders(1:end-1), bhv.movingFracNoGo(noGoDay(iDay), :), 'k');
        p3.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
        ylim([0 1])
    end
    ylim([0 1])
    xlim([binBorders(1), binBorders(end)])
    legend([p3], {'NoGo'})
    xlabel('time from stim onset (s)')
    ylabel('fraction moving')
    makepretty;
    ylim([0 1])
    end
    
    %an(iAnimal).movingFracGo = movingFracGo;
    %an(iAnimal).movingFracGo = movingFracNoGo;
    an(iAnimal).noGo = bhv.noGo;
    an(iAnimal).nTrials = bhv.nTrials;
    an(iAnimal).goLeft = bhv.goLeft;
    an(iAnimal).noGoDay = noGoDay;
end

% figure();
% plot(binBorders(1:end-1), an(1).movingFrac);
% hold on;
% makepretty;
% plot(binBorders(1:end-1), an(2).movingFrac);
% hold on;
% makepretty;
% plot(binBorders(1:end-1), an(3).movingFrac);
% hold on;
% makepretty;
% plot(binBorders(1:end-1), an(4).movingFrac);
% hold on;
% legend(animalsAll)
% xlabel('time from stim onset (s)')
% ylabel('fraction moving')
% makepretty;

%% summary: % no go per stim per mouse on last day only 

figure();

for iMouse =[5:6]
subplot(2,1,iMouse-4)
transparencyValues = 0:1 / length(an(iMouse).noGoDay):1;
for iDay = 1:size(an(iMouse).noGoDay,2)
    if iDay == size(an(iMouse).noGoDay,2)
        col1 = rgb('DarkBlue');
        col2 = rgb('DarkRed');
        
    else
        col1 = 'b';
        col2='r';
    end
    scatter([1,2,3],an(iMouse).goLeft(an(iMouse).noGoDay(iDay), [1,3,2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1,3,2]),[],col1,'filled',...
        'MarkerFaceAlpha',transparencyValues(iDay+1),'MarkerEdgeAlpha',transparencyValues(iDay+1))
    p=plot([1,2,3],an(iMouse).goLeft(an(iMouse).noGoDay(iDay), [1,3,2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1,3,2]),'Color', col1);
    p.Color(4)=transparencyValues(iDay+1);
    hold on;
    makepretty;
    scatter([1,2,3],an(iMouse).noGo(an(iMouse).noGoDay(iDay), [1,3,2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1,3,2]),[],col2,'filled',...
        'MarkerFaceAlpha',transparencyValues(iDay+1),'MarkerEdgeAlpha',transparencyValues(iDay+1))
     p=plot([1,2,3],an(iMouse).noGo(an(iMouse).noGoDay(iDay), [1,3,2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1,3,2]),'Color', col2);
    p.Color(4)=transparencyValues(iDay+1);
    makepretty;
end

xlabel('image type')
ylabel('frac go left(blue)/nogo (red)')
xticks([1 2 3])
xticklabels({'Go1','NoGo','Go2'})
ylim([0 1])
makepretty;
grid on;
end


