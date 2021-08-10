clear all;
%clc;

plotMore=false;
% to add first half, second half
% add weight %
% 2021 04 08 added no go analysis
animalsAll = {'JF036', 'JF037', 'JF038', 'JF039', 'JF042', 'JF043', 'JF044'};
animalsPhase1 = {[1, 2], [1, 2], [1:4]};
animalsPhase2 = {[3, 4], [3:5], [5:9]};
animalsPhase3 = {[5], [6], [10]};
animalsPhase4 = {[6:12], [7:14], [11:13]};
animalsPhase5 = {[13, 14, 16], [16], []};

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
    for curr_day = theseD %[1:19, 21:length(experiments)]%JF039


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
                (block.events.hitValues(response_trials) == 0) & go1Trials ...
                );
            noGo((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                (block.events.hitValues(response_trials) == 0) & go2Trials ...
                );
            noGo((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                (block.events.hitValues(response_trials) == 0) & noGoTrials ...
                );

            noGo1((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials(1:round(size(response_trials, 2)/2))) & ...
                block.events.trialSideValues(response_trials(1:round(size(response_trials, 2)/2))) == 1 & ...
                (block.events.hitValues(response_trials(1:round(size(response_trials, 2)/2))) == 0) & go1Trials(1:round(size(response_trials, 2)/2)) ...
                );
            noGo1((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials(1:round(size(response_trials, 2)/2))) & ...
                block.events.trialSideValues(response_trials(1:round(size(response_trials, 2)/2))) == 1 & ...
                (block.events.hitValues(response_trials(1:round(size(response_trials, 2)/2))) == 0) & go2Trials(1:round(size(response_trials, 2)/2)) ...
                );
            noGo1((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials(1:round(size(response_trials, 2)/2))) & ...
                block.events.trialSideValues(response_trials(1:round(size(response_trials, 2)/2))) == 1 & ...
                (block.events.hitValues(response_trials(1:round(size(response_trials, 2)/2))) == 0) & noGoTrials(1:round(size(response_trials, 2)/2)) ...
                );

            noGo2((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials(round(size(response_trials, 2)/2):end)) & ...
                block.events.trialSideValues(response_trials(round(size(response_trials, 2)/2):end)) == 1 & ...
                (block.events.hitValues(response_trials(round(size(response_trials, 2)/2):end)) == 0) & go1Trials(round(size(response_trials, 2)/2):end) ...
                );
            noGo2((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials(round(size(response_trials, 2)/2):end)) & ...
                block.events.trialSideValues(response_trials(round(size(response_trials, 2)/2):end)) == 1 & ...
                (block.events.hitValues(response_trials(round(size(response_trials, 2)/2):end)) == 0) & go2Trials(round(size(response_trials, 2)/2):end) ...
                );
            noGo2((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials(round(size(response_trials, 2)/2):end)) & ...
                block.events.trialSideValues(response_trials(round(size(response_trials, 2)/2):end)) == 1 & ...
                (block.events.hitValues(response_trials(round(size(response_trials, 2)/2):end)) == 0) & noGoTrials(round(size(response_trials, 2)/2):end) ...
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

            %trial_stim = block.events.trialContrastValues(response_trials).*block.events.trialSideValues(response_trials);
            %stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1;1],[],1));
            %[~,trial_stim_idx] = ismember(trial_stim,stim_list);
            if isfield(block.events, 'stimOnTimes')
                trial_wheel_starts = arrayfun(@(x) ...
                    wheel_starts(find(wheel_starts > block.events.stimOnTimes(x), 1)), ...
                    response_trials);

                trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
                stim_rxn_time = nanmedian(trial_move_t);
                stim_rxn_timeSEM = nanstd(trial_move_t) ./ sqrt(numel(trial_move_t));
            else
                try
                    trial_wheel_starts = arrayfun(@(x) ...
                        wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                        response_trials(1:end));
                    trial_wheel_move_start = arrayfun(@(x) ...
                        find(wheel_t_resample > block.events.stimulusOnTimes(x), 1), ...
                        response_trials(1:end));
                    trial_wheel_move = nan(size(response_trials(1:end),2), 1600); %about 1.6 seconds 
                    for iTrial = 1:size(response_trials(1:end),2)
                        trial_wheel_move(iTrial, :) = wheel_velocity(trial_wheel_move_start(iTrial):...
                            trial_wheel_move_start(iTrial)+1599);
                    end
                    
                    trial_move_t = trial_wheel_starts - block.events.stimulusOnTimes(response_trials(1:end));
                    stim_rxn_time = nanmedian(trial_move_t);
                    stim_rxn_timeSEM = nanstd(trial_move_t) ./ sqrt(numel(trial_move_t));
                catch
                    trial_wheel_starts = NaN;
                    trial_move_t = NaN;
                    stim_rxn_time = NaN;
                    stim_rxn_timeSEM = NaN;
                    trial_wheel_move = nan(size(response_trials(1:end),2), 1600);
                end
            end

            %useless
            %             forWS = [0, block.events.stimOffTimes];
            %             for iStim = 1:length(block.events.stimOffTimes)
            %                 wheel_startsITI(iStim) = numel(sum(wheel_starts > forWS(iStim) &  wheel_starts< block.events.stimOnTimes(iStim)));
            %                 wheel_startsStim(iStim) = numel(sum(wheel_starts > block.events.stimOnTimes(iStim) &  wheel_starts< block.events.stimOffTimes(iStim)));
            %             end

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


                if sum(go2Trials) >= 1
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
                    movingFracNoGo = nan(1, 250);
                end
            end


            bhv.correctTrials(curr_day, :) = correctTrials;
            bhv.nTrials(curr_day, :) = nTrials;
            bhv.n_trials(curr_day, :) = n_trials;
            bhv.total_water(curr_day, :) = total_water;
            bhv.conditions(curr_day, :) = conditi;
            bhv.goLeft(curr_day, :) = goLeft;
            bhv.noGo(curr_day, :) = noGo;
            bhv.noGo1(curr_day, :) = noGo1;
            bhv.noGo2(curr_day, :) = noGo2;
            bhv.trialTime(curr_day, :) = trialTime;
            
            keep_day = [keep_day, curr_day];
            if contains(block.expDef, 'NoGo')
                noGoDay = [noGoDay, curr_day];
                if sum(go2Trials) >= 1

                    bhv.movingFracGo2(curr_day, :) = movingFracGo2;
                end
                bhv.movingFracGo1(curr_day, :) = movingFracGo1;
                bhv.movingFracNoGo(curr_day, :) = movingFracNoGo;
                bhv.trialMoveNoGoCorrect{curr_day,:,:} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0),:);
                bhv.trialMoveNoGoIncorrect{curr_day,:,:} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1),:);
                bhv.trialMoveGo1Correct{curr_day,:,:} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1),:);
                bhv.trialMoveGo2Correct{curr_day,:,:} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1),:);
                 bhv.trialMoveGo1Incorrect{curr_day,:,:} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0),:);
                bhv.trialMoveGo2Incorrect{curr_day,:,:} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0),:);
            else
                bhv.movingFrac(curr_day, :) = movingFrac;
              
               

            end
            bhv.stim_rxn_time(curr_day) = stim_rxn_time;
            bhv.stim_rxn_timeSEM(curr_day) = stim_rxn_timeSEM;

            %stims
        end
    end

    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels = cellfun(@(day) [' ', day(6:end)], ...
        {experiments(keep_day).day}, 'uni', false);
    %
    figure('Name', animal);
    subplot(331)
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
    % figure();
    subplot(333)
    con = bhv.conditions(keep_day(1), :);
    ccc = num2str(con(1, :));
    cc = textscan(ccc, '%s', 'Delimiter', ' ');
    cc = cc{1};
    cc(strcmp('', cc)) = [];
    im = imagesc(1:length(con), 1:size(bhv.correctTrials(keep_day, :), 1), bhv.goLeft(keep_day, :)./bhv.nTrials(keep_day, :));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], 'RdBu'));
    c = colorbar;
    ylabel(c, 'Go left (frac)');
    xlabel('image type');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);
    xticks([1, 2, 3])
    xticklabels({'Go1', 'Go2', 'NoGo'})
    %
    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.correctTrials(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    %axis square;
    caxis([0, 1])
    makepretty;
    %  figure();
    subplot(334)
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
    xticks([1, 2, 3])
    xticklabels({'Go1', 'Go2', 'NoGo'})

    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.noGo(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    %axis square;
    caxis([0, 1])
    makepretty;

    subplot(332)
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

    subplot(335)

    transparencyValues = 0:1 / length(noGoDay):1;



    if length(noGoDay) == 0
        transparencyValues = 0:1 / length(keep_day):1;
        for iDay = 1:length(keep_day)
            p1 = plot(binBorders(1:end-1), bhv.movingFrac(keep_day(iDay), :), 'b');
            p1.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
        end
        ylim([0, 1])
    else
        for iDay = 1:length(noGoDay)
            p1 = plot(binBorders(1:end-1), bhv.movingFracGo1(noGoDay(iDay), :), 'b');
            p1.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
            ylim([0, 1])

        end
        if isfield(bhv, 'movingFracGo2')
            p2 = plot(binBorders(1:end-1), bhv.movingFracGo2(noGoDay(iDay), :), 'r');
            p2.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
            legend([p1, p2], {'Go1', 'Go2'})
        end
        ylim([0, 1])
    end
    xlim([binBorders(1), binBorders(end)])

    xlabel('time from stim onset (s)')
    ylabel('fraction moving')
    makepretty;

    subplot(336)
    if length(noGoDay) == 0
    else
        transparencyValues = 0:1 / length(noGoDay):1;
        for iDay = 1:length(noGoDay)

            p3 = plot(binBorders(1:end-1), bhv.movingFracNoGo(noGoDay(iDay), :), 'k');
            p3.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
            ylim([0, 1])
        end
        ylim([0, 1])
        xlim([binBorders(1), binBorders(end)])
        legend([p3], {'NoGo'})
        xlabel('time from stim onset (s)')
        ylabel('fraction moving')
        makepretty;
        ylim([0, 1])
    end
    for iNoGoday = 1:size(noGoDay,2)
        meanGo1Correct(iNoGoday,:) = nanmean(abs(bhv.trialMoveGo1Correct{noGoDay(iNoGoday),:,:}));
        meanGo2Correct(iNoGoday,:) = nanmean(abs(bhv.trialMoveGo2Correct{noGoDay(iNoGoday),:,:}));
        meanNoGoCorrect(iNoGoday,:) = nanmean(abs(bhv.trialMoveNoGoCorrect{noGoDay(iNoGoday),:,:}));
        meanGo1Incorrect(iNoGoday,:) = nanmean(abs(bhv.trialMoveGo1Incorrect{noGoDay(iNoGoday),:,:}));
        meanGo2Incorrect(iNoGoday,:) = nanmean(abs(bhv.trialMoveGo2Incorrect{noGoDay(iNoGoday),:,:}));
        meanNoGoIncorrect(iNoGoday,:) = nanmean(abs(bhv.trialMoveNoGoIncorrect{noGoDay(iNoGoday),:,:}));
    end
    subplot(337) 
    cla;
    plot(0:1.6/1600:1.6-1.6/1600,nanmean( meanGo1Correct),'Color', 'b'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd( meanGo1Correct) + ...
        nanmean( meanGo1Correct); nanstd( meanGo1Correct) + nanmean( meanGo1Correct)],'b')
    makepretty; 
    
   
    plot(0:1.6/1600:1.6-1.6/1600,nanmean( meanGo2Correct),'Color', 'r'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd( meanGo2Correct) + ...
        nanmean( meanGo2Correct); nanstd( meanGo2Correct) + nanmean( meanGo2Correct)],'r')
    makepretty; 
    
  
    plot(0:1.6/1600:1.6-1.6/1600,nanmean(meanNoGoIncorrect),'Color', 'k'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd(meanNoGoIncorrect) + ...
        nanmean(meanNoGoIncorrect); nanstd(meanNoGoIncorrect) + nanmean(meanNoGoIncorrect)],'k')
    makepretty; 
    xlim([0 1.6])
    
    subplot(338)
    cla;
    plot(0:1.6/1600:1.6-1.6/1600,nanmean( meanGo1Incorrect),'Color', 'b'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd( meanGo1Incorrect) + ...
        nanmean( meanGo1Incorrect); nanstd( meanGo1Incorrect) + nanmean( meanGo1Incorrect)],'b')
    makepretty; 

       plot(0:1.6/1600:1.6-1.6/1600,nanmean( meanGo2Incorrect),'Color', 'r'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd( meanGo2Incorrect) + ...
        nanmean( meanGo2Incorrect); nanstd( meanGo2Incorrect) + nanmean( meanGo2Incorrect)],'r')
    makepretty; 
    
  
    plot(0:1.6/1600:1.6-1.6/1600,nanmean(meanNoGoCorrect),'Color', 'k'); hold on;
    plotshaded(0:1.6/1600:1.6-1.6/1600, [- nanstd(meanNoGoCorrect) + ...
        nanmean(meanNoGoCorrect); nanstd(meanNoGoCorrect) + nanmean(meanNoGoCorrect)],'k')
    makepretty; 
    xlim([0 1.6])

    
    
%     subplot(326)
%     if length(noGoDay)==0
%         else
%     transparencyValues = 0:1 / length(noGoDay):1;
%     for iDay = 1:length(noGoDay)
%         
%         p3 = plot(binBorders(1:end-1), bhv.movingFracNoGo(noGoDay(iDay), :), 'k');
%         p3.Color(4) = transparencyValues(iDay+1);
%         makepretty;
%         hold on;
%         ylim([0 1])
%     end
%     ylim([0 1])
%     xlim([binBorders(1), binBorders(end)])
%     legend([p3], {'NoGo'})
%     xlabel('time from stim onset (s)')
%     ylabel('fraction moving')
%     makepretty;
%     ylim([0 1])
%     end

    

    
    %plot(nanmean(bhv.
    %an(iAnimal).movingFracGo = movingFracGo;
    %an(iAnimal).movingFracGo = movingFracNoGo;
    an(iAnimal).noGo = bhv.noGo;
    an(iAnimal).noGo1 = bhv.noGo1;
    an(iAnimal).noGo2 = bhv.noGo2;
    an(iAnimal).nTrials = bhv.nTrials;
    an(iAnimal).goLeft = bhv.goLeft;
    an(iAnimal).noGoDay = noGoDay;

    an(iAnimal).stim_rxn_time = bhv.stim_rxn_time;
    an(iAnimal).movingFrac = bhv.movingFrac;
    an(iAnimal).n_trials_all = bhv.n_trials;
    clearvars bhv

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

for iAnimal = [5:7]
    subplot(3, 1, iAnimal-4)
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    for iDay = 1:size(an(iAnimal).noGoDay, 2)
        if iDay == size(an(iAnimal).noGoDay, 2)
            col1 = rgb('DarkBlue');
            col2 = rgb('DarkRed');

        else
            col1 = 'b';
            col2 = 'r';
        end
        scatter([1, 2, 3], an(iAnimal).goLeft(an(iAnimal).noGoDay(iDay), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 3, 2]), [], col1, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).goLeft(an(iAnimal).noGoDay(iDay), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 3, 2]), 'Color', col1);
        p.Color(4) = transparencyValues(iDay+1);
        hold on;
        makepretty;
        scatter([1, 2, 3], an(iAnimal).noGo(an(iAnimal).noGoDay(iDay), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 3, 2]), [], col2, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).noGo(an(iAnimal).noGoDay(iDay), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 3, 2]), 'Color', col2);
        p.Color(4) = transparencyValues(iDay+1);
        makepretty;
    end

    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')

    xticks([1, 2, 3])
    xticklabels({'Go1', 'NoGo', 'Go2'})
    ylim([0, 1])
    makepretty;
    grid on;
end

figure();

for iAnimal = [5:7]
    subplot(3, 1, iAnimal-4)
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    for iDay = size(an(iAnimal).noGoDay, 2)
        if iDay == size(an(iAnimal).noGoDay, 2)
            col1 = rgb('DarkBlue');
            col2 = rgb('DarkRed');

        else
            col1 = 'b';
            col2 = 'r';
        end
        scatter([1, 2, 3], an(iAnimal).goLeft(an(iAnimal).noGoDay(iDay), [1, 2, 3])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 2, 3]), [], col1, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).goLeft(an(iAnimal).noGoDay(iDay), [1, 2, 3])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 2, 3]), 'Color', col1);
        p.Color(4) = transparencyValues(iDay+1);
        hold on;
        makepretty;
        scatter([1, 2, 3], an(iAnimal).noGo(an(iAnimal).noGoDay(iDay), [1, 2, 3])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 2, 3]), [], col2, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).noGo(an(iAnimal).noGoDay(iDay), [1, 2, 3])./an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), [1, 2, 3]), 'Color', col2);
        p.Color(4) = transparencyValues(iDay+1);
        makepretty;
    end

    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')

    xticks([1, 2, 3])
    xticklabels({'Go1', 'Go2', 'NoGo'})
    ylim([0, 1])
    makepretty;
    grid on;
end





if plotMore == true
    figure();

for iAnimal = 5:6
    subplot(2, 1, iAnimal-4)
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    if iAnimal == 6
        thiD = size(an(iAnimal).noGoDay, 2) -2:size(an(iAnimal).noGoDay, 2);
        % thiD=thiD(1,2,4);
    elseif iAnimal == 5
        thiD = size(an(iAnimal).noGoDay, 2) - 2:size(an(iAnimal).noGoDay, 2);
        %thiD=thiD([1,2,4,5]);
    else


    end


    col1 = 'b';
    col2 = 'r';
    if length(thiD) == 1
        thisS = [an(iAnimal).noGo(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))];
    else
        thisS = nanmean([an(iAnimal).noGo(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
        thisSt = nanstd([an(iAnimal).noGo(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    end
    scatter([1, 2, 3], thisS, [], col2, 'filled')
    hold on;
    if length(thiD) == 1
    else
        errorbar([1, 2, 3], thisS, thisSt, thisSt, 'Color', col2, 'LineWidth', 2)
    end
    p1 = plot([1, 2, 3], thisS, 'Color', col2);
    makepretty;
    hold on;
    if length(thiD) == 1
        thisS = [an(iAnimal).goLeft(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))];
    else

        thisS = nanmean([an(iAnimal).goLeft(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
        thisSt = nanstd([an(iAnimal).goLeft(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    end
    scatter([1, 2, 3], thisS, [], col1, 'filled')
    hold on;
    if length(thiD) == 1
    else
        errorbar([1, 2, 3], thisS, thisSt, thisSt, 'Color', col1, 'LineWidth', 2)
    end

    p2 = plot([1, 2, 3], thisS, 'Color', col1);
    makepretty;
    hold on;

    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
    xticks([1, 2, 3])
    xticklabels({'Go1', 'Go2', 'NoGo'})
    ylim([0, 1])
    makepretty;
    grid on;


end
xlabel('image type')
ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
xticks([1, 2, 3])
xticklabels({'Go1', 'Go2', 'NoGo'})
ylim([0, 1])
makepretty;
grid on;

figure();

for iAnimal = [5:7]
    subplot(3, 2, (iAnimal - 4 - 1)*2+1)
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    for iDay = 1:size(an(iAnimal).noGoDay, 2)
        if iDay == size(an(iAnimal).noGoDay, 2)
            col1 = rgb('DarkBlue');
            col2 = rgb('DarkRed');

        else
            col1 = 'b';
            col2 = 'r';
        end
        thisS = [an(iAnimal).noGo1(an(iAnimal).noGoDay(iDay), 1) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), 1) ./ 2), ...
            an(iAnimal).noGo2(an(iAnimal).noGoDay(iDay), 1) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), 1) ./ 2)];
        scatter([1, 2], thisS, [], col2, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2], thisS, 'Color', col2);
        p.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
    end


    xlabel('image type')
    ylabel('frac nogo on go1 trials(red)')
    xticks([1, 2, 3])
    xticklabels({'Go1 1/2 trials', 'Go1 2/2 trials'})
    ylim([0, 1])
    makepretty;
    grid on;
    subplot(3, 2, (iAnimal - 4 - 1)*2+2)
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    for iDay = 1:size(an(iAnimal).noGoDay, 2)
        if iDay == size(an(iAnimal).noGoDay, 2)
            col1 = rgb('DarkBlue');
            col2 = rgb('DarkRed');

        else
            col1 = 'b';
            col2 = 'r';
        end
        thisS = [an(iAnimal).noGo1(an(iAnimal).noGoDay(iDay), 3) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), 3) ./ 2), ...
            an(iAnimal).noGo2(an(iAnimal).noGoDay(iDay), 3) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(iDay), 3) ./ 2)];
        scatter([1, 2], thisS, [], col2, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2], thisS, 'Color', col2);
        p.Color(4) = transparencyValues(iDay+1);
        makepretty;
        hold on;
    end


    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')

    xticks([1, 2, 3])
    xticklabels({'NoGo 1/2 trials', 'NoGo 2/2 trials'})
    ylim([0, 1])
    makepretty;
    grid on;

end
% plot(1:10, rand(1,10))
% ax = gca;
% % Simply color an XTickLabel
% ax.XTickLabel{3} = ['\color{red}' ax.XTickLabel{3}];
% % Use TeX symbols
% ax.XTickLabel{4} = '\color{blue} \uparrow';
% % Use multiple colors in one XTickLabel
% ax.XTickLabel{5} = '\color[rgb]{0,1,0}green\color{orange}?';
% % Color YTickLabels with colormap
% nColors = numel(ax.YTickLabel);
% cm = jet(nColors);
% for i = 1:nColors
%     ax.YTickLabel{i} = sprintf('\color[rgb]{%f,%f,%f}%s', cm(i,:), ax.YTickLabel{i});
% end

%% Going Left phase 1
figure();
animalCol = {rgb('HotPink'); rgb('Crimson'); rgb('RoyalBlue')};
for iAnimal = 5:7
    cmap = brewermap(length(animalsPhase1{iAnimal-4})+2, 'Greens');

    scatter(1:length(animalsPhase1{iAnimal-4}), ...
        an(iAnimal).goLeft(animalsPhase1{iAnimal-4}, 1)./an(iAnimal).nTrials(animalsPhase1{iAnimal-4}, 1), [], animalCol{iAnimal-4}, 'filled');
    hold on;
    p(iAnimal) = plot(1:length(animalsPhase1{iAnimal-4}), ...
        an(iAnimal).goLeft(animalsPhase1{iAnimal-4}, 1)./an(iAnimal).nTrials(animalsPhase1{iAnimal-4}, 1), 'Color', animalCol{iAnimal-4});
    makepretty;
    hold on;

end
legend([p(5), p(6), p(7)], {'JF042', 'JF043', 'JF044'})
ylabel('Frac. go left')
xticks([1, 2, 3, 4])
xlabel('Day #')
grid on;
makeprettyLarge;

%% Fraction moving phase 1 vs 2
figure();

for iAnimal = 5:7
    meanMovingFrac(iAnimal, :) = nanmean(an(iAnimal).movingFrac(animalsPhase1{iAnimal-4}, :));
    %     plot( meanMovingFrac(iAnimal,:), 'Color',rgb('DeepSkyBlue'))
    %     makepretty;
    %      hold on;
end
plot(binBorders(1:end-1), nanmean(meanMovingFrac), 'Color', rgb('RoyalBlue'))
fuckyou = nanmean(meanMovingFrac) - (nanstd(meanMovingFrac) ./ sqrt(3));
fuckyoutoo = nanmean(meanMovingFrac) + (nanstd(meanMovingFrac) ./ sqrt(3));
plotshaded(binBorders(1:end-1), [fuckyou; fuckyoutoo], rgb('RoyalBlue'));
makepretty;
hold on;
for iAnimal = 5:7
    meanMovingFrac(iAnimal, :) = nanmean(an(iAnimal).movingFrac(animalsPhase2{iAnimal-4}(1:end - 1), :));
    %     plot( meanMovingFrac(iAnimal,:), 'Color',rgb('DeepSkyBlue'))
    %     makepretty;
    %      hold on;
end
plot(binBorders(1:end-1), nanmean(meanMovingFrac), 'Color', rgb('Green'))
fuckyou = nanmean(meanMovingFrac) - (nanstd(meanMovingFrac) ./ sqrt(3));
fuckyoutoo = nanmean(meanMovingFrac) + (nanstd(meanMovingFrac) ./ sqrt(3));
plotshaded(binBorders(1:end-1), [fuckyou; fuckyoutoo], rgb('Green'));
makepretty;
hold on;
for iAnimal = 5:7
    meanMovingFrac(iAnimal, :) = an(iAnimal).movingFrac(animalsPhase2{iAnimal-4}(end), :);
    %     plot( meanMovingFrac(iAnimal,:), 'Color',rgb('DeepSkyBlue'))
    %     makepretty;
    %      hold on;
end
plot(binBorders(1:end-1), nanmean(meanMovingFrac), 'Color', rgb('Purple'))
fuckyou = nanmean(meanMovingFrac) - (nanstd(meanMovingFrac) ./ sqrt(3));
fuckyoutoo = nanmean(meanMovingFrac) + (nanstd(meanMovingFrac) ./ sqrt(3));
plotshaded(binBorders(1:end-1), [fuckyou; fuckyoutoo], rgb('Purple'));
makepretty;
xlim([binBorders(1), binBorders(end)])
xlabel('time from stim onset (s)')
ylabel('fraction moving')
L(1) = plot(nan, nan, 'Color', rgb('RoyalBlue'), 'LineWidth', 2);
L(2) = plot(nan, nan, 'Color', rgb('Green'), 'LineWidth', 2);
L(3) = plot(nan, nan, 'Color', rgb('Purple'), 'LineWidth', 2);
legend(L, {'Phase 1', 'Phase 2', 'Last day phase 2'})

makeprettyLarge;

%% moving frac phase 2 vs 3
figure();
for iAnimal = 5:7
    meanMovingFrac(iAnimal, :) = an(iAnimal).movingFrac(animalsPhase3{iAnimal-4}(end), :);
    %     plot( meanMovingFrac(iAnimal,:), 'Color',rgb('DeepSkyBlue'))
    %     makepretty;
    %      hold on;
end
plot(binBorders(1:end-1), nanmean(meanMovingFrac), 'Color', rgb('Orange'))
fuckyou = nanmean(meanMovingFrac) - (nanstd(meanMovingFrac) ./ sqrt(3));
fuckyoutoo = nanmean(meanMovingFrac) + (nanstd(meanMovingFrac) ./ sqrt(3));
plotshaded(binBorders(1:end-1), [fuckyou; fuckyoutoo], rgb('Orange'));
makepretty;
hold on;
for iAnimal = 5:7
    meanMovingFrac(iAnimal, :) = an(iAnimal).movingFrac(animalsPhase2{iAnimal-4}(end), :);
    %     plot( meanMovingFrac(iAnimal,:), 'Color',rgb('DeepSkyBlue'))
    %     makepretty;
    %      hold on;
end
plot(binBorders(1:end-1), nanmean(meanMovingFrac), 'Color', rgb('Purple'))
fuckyou = nanmean(meanMovingFrac) - (nanstd(meanMovingFrac) ./ sqrt(3));
fuckyoutoo = nanmean(meanMovingFrac) + (nanstd(meanMovingFrac) ./ sqrt(3));
plotshaded(binBorders(1:end-1), [fuckyou; fuckyoutoo], rgb('Purple'));
makepretty;
xlim([binBorders(1), binBorders(end)])
xlabel('time from stim onset (s)')
ylabel('fraction moving')
L(1) = plot(nan, nan, 'Color', rgb('Purple'), 'LineWidth', 2);
L(2) = plot(nan, nan, 'Color', rgb('Orange'), 'LineWidth', 2);
legend(L, {'Last day phase 2', 'Phase 3'})
makeprettyLarge;

%% Reaction times phase 1 vs 2
figure();
for iAnimal = 5:7
    cmap = brewermap(length(animalsPhase1{iAnimal-4})+3, 'Greens');
    errorbar(iAnimal-4, nanmean(an(iAnimal).stim_rxn_time(animalsPhase1{iAnimal-4})), ...
        nanstd(an(iAnimal).stim_rxn_time(animalsPhase1{iAnimal-4})/sqrt(length(an(iAnimal).stim_rxn_time(animalsPhase1{iAnimal-4})))), ...
        'Color', 'k', 'LineWidth', 2);
    makepretty;
    hold on;
    scatter(repmat(iAnimal-4, length(animalsPhase1{iAnimal-4}), 1), an(iAnimal).stim_rxn_time(animalsPhase1{iAnimal-4}), [], cmap(2:end-2, :), 'filled');
    makepretty;
    hold on;
    scatter(iAnimal-4, nanmean(an(iAnimal).stim_rxn_time(animalsPhase1{iAnimal-4})), [], 'k', 'filled');
    makepretty;
    hold on;
end
ylabel('Reaction time (s)')
xticks([1, 2, 3])
xticklabels({'JF042', 'JF043', 'JF044'})
legend({'mean +/ - se'})
ylim([0, 6])
grid on;


for iAnimal = 5:7
    cmap = brewermap(length(animalsPhase2{iAnimal-4})+3, 'Blues');

    errorbar(iAnimal-4+0.1, nanmean(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})), ...
        nanstd(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})/sqrt(length(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})))), ...
        'Color', 'k', 'LineWidth', 2);
    makepretty;
    hold on;
    scatter(repmat(iAnimal-4+0.1, length(animalsPhase2{iAnimal-4}), 1), an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4}), [], cmap(2:end-2, :), 'filled');
    makepretty;
    hold on;
    scatter(iAnimal-4+0.1, nanmean(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})), [], 'k', 'filled');
    makepretty;
    hold on;
end
ylabel('Reaction time (s)')
xticks([1, 2, 3])
xticklabels({'JF042', 'JF043', 'JF044'})
L(3) = scatter(nan, nan, [], rgb('DeepSkyBlue'), 'filled');
L(2) = scatter(nan, nan, [], rgb('MediumSeaGreen'), 'filled');
L(1) = errorbar(nan, nan, 'Color', 'k', 'LineWidth', 2);
legend(L, {'mean +/ - SE', 'Phase 1', 'Phase 2'})
%legend({'mean +/ - SE'})
ylim([0, 6])
grid on;
xlim([0.5, 3.6])
makeprettyLarge;

%% reaction time phase 2 and 3
figure();
for iAnimal = 5:7
    cmap = brewermap(length(animalsPhase2{iAnimal-4})+3, 'Blues');

    errorbar(iAnimal-4, nanmean(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})), ...
        nanstd(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})/sqrt(length(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})))), ...
        'Color', 'k', 'LineWidth', 2);
    makepretty;
    hold on;
    scatter(repmat(iAnimal-4, length(animalsPhase2{iAnimal-4}), 1), an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4}), [], cmap(2:end-2, :), 'filled');
    makepretty;
    hold on;
    scatter(iAnimal-4, nanmean(an(iAnimal).stim_rxn_time(animalsPhase2{iAnimal-4})), [], 'k', 'filled');
    makepretty;
    hold on;
end
ylabel('Reaction time (s)')
xticks([1, 2, 3])
xticklabels({'JF042', 'JF043', 'JF044'})
legend({'mean +/ - se'})
ylim([0, 6])
grid on;


for iAnimal = 5:7
    cmap = brewermap(length(animalsPhase3{iAnimal-4})+3, 'Oranges');

    errorbar(iAnimal-4+0.1, nanmean(an(iAnimal).stim_rxn_time(animalsPhase3{iAnimal-4})), ...
        nanstd(an(iAnimal).stim_rxn_time(animalsPhase3{iAnimal-4})/sqrt(length(an(iAnimal).stim_rxn_time(animalsPhase3{iAnimal-4})))), ...
        'Color', 'k', 'LineWidth', 2);
    makepretty;
    hold on;
    scatter(repmat(iAnimal-4+0.1, length(animalsPhase3{iAnimal-4}), 1), an(iAnimal).stim_rxn_time(animalsPhase3{iAnimal-4}), [], cmap(3:end-1, :), 'filled');
    makepretty;
    hold on;
    %scatter(iAnimal-4+0.1, nanmean(an(iAnimal).stim_rxn_time(animalsPhase3{iAnimal-4})),[],'k', 'filled') ;
    makepretty;
    hold on;
end
ylabel('Reaction time (s)')
xticks([1, 2, 3])
xticklabels({'JF042', 'JF043', 'JF044'})
L(2) = scatter(nan, nan, [], rgb('DeepSkyBlue'), 'filled');
L(3) = scatter(nan, nan, [], cmap(3:end-1, :), 'filled');
L(1) = errorbar(nan, nan, 'Color', 'k', 'LineWidth', 2);
legend(L, {'mean +/ - SE', 'Phase 2', 'Phase 3'})
%legend({'mean +/ - SE'})
ylim([0, 2.1])
grid on;
xlim([0.5, 3.6])
makeprettyLarge;

%% Phase 4: no going

figure();
for iAnimal = 5:7
    for iDay = 1:length(animalsPhase4{iAnimal-4})
        cmap = copper(length(animalsPhase4{iAnimal-4}+1));
        plot([iAnimal - 4, iAnimal - 4 + 0.5], [an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(iDay)) ./ an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay)), ...
            an(iAnimal).noGo(animalsPhase4{iAnimal-4}(iDay)) ./ an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay))], 'Color', cmap(size(cmap, 1)+1-iDay, :));
        makepretty;
        hold on;
        scatter([iAnimal - 4, iAnimal - 4 + 0.5], [an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(iDay)) ./ an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay)), ...
            an(iAnimal).noGo(animalsPhase4{iAnimal-4}(iDay)) ./ an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay))], [], cmap(size(cmap, 1)+1-iDay, :), 'filled')
    end
end

figure();

for iAnimal = 5
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    if iAnimal == 7

        thiD = size(an(iAnimal).noGoDay, 2) - 1:size(an(iAnimal).noGoDay, 2);
    elseif iAnimal == 5
        thiD = size(an(iAnimal).noGoDay, 2) - 1:size(an(iAnimal).noGoDay, 2);
    else


    end


    col1 = 'b';
    col2 = 'r';

    thisS = nanmean([an(iAnimal).noGo(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    thisSt = nanstd([an(iAnimal).noGo(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    scatter([1, 2, 3], thisS, [], col2, 'filled')
    hold on;
    errorbar([1, 2, 3], thisS, thisSt, thisSt, 'Color', col2, 'LineWidth', 2)
    p1 = plot([1, 2, 3], thisS, 'Color', col2);
    makepretty;
    hold on;

    thisS = nanmean([an(iAnimal).goLeft(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    thisSt = nanstd([an(iAnimal).goLeft(an(iAnimal).noGoDay(thiD), [1, 2, 3]) ./ (an(iAnimal).nTrials(an(iAnimal).noGoDay(thiD), [1, 2, 3]))]);
    scatter([1, 2, 3], thisS, [], col1, 'filled')
    hold on;
    errorbar([1, 2, 3], thisS, thisSt, thisSt, 'Color', col1, 'LineWidth', 2)
    p2 = plot([1, 2, 3], thisS, 'Color', col1);
    makepretty;
    hold on;


end
xlabel('image type')
ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
xticks([1, 2, 3])
xticklabels({'Go1', 'Go2', 'NoGo'})
ylim([0, 1])
makepretty;
grid on;

figure();

for iAnimal = [5:7]
    subplot(3, 1, iAnimal-4)
    transparencyValues = 0:1 / size(animalsPhase4{iAnimal-4}, 2):1;
    for iDay = 1:size(animalsPhase4{iAnimal-4}, 2)
        if iDay == size(animalsPhase4{iAnimal-4}, 2)
            col1 = rgb('DarkBlue');
            col2 = rgb('DarkRed');

        else
            col1 = 'b';
            col2 = 'r';
        end
        scatter([1, 2, 3], an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2])./an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2]), [], col1, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2])./an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2]), 'Color', col1);
        p.Color(4) = transparencyValues(iDay+1);
        hold on;
        makepretty;
        scatter([1, 2, 3], an(iAnimal).noGo(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2])./an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2]), [], col2, 'filled', ...
            'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
        p = plot([1, 2, 3], an(iAnimal).noGo(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2])./an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(iDay), [1, 3, 2]), 'Color', col2);
        p.Color(4) = transparencyValues(iDay+1);
        makepretty;
    end

    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')

    xticks([1, 2, 3])
    xticklabels({'Go1', 'NoGo', 'Go2'})
    ylim([0, 1])
    makepretty;
    grid on;
end

figure(); %last two days
iAnimal = 6
plot(nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), 'Color', rgb('RoyalBlue'))
makepretty;
hold on;
errorbar([1, 2], nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), ...
    nanstd([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]./sqrt(2)), ...
    'Color', rgb('RoyalBlue'), 'LineWidth', 2)
makepretty;
hold on;
plot(nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), 'Color', rgb('IndianRed'))
makepretty;
hold on;
errorbar([1, 2], nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), ...
    nanstd([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]./sqrt(2)), ...
    'Color', rgb('IndianRed'), 'LineWidth', 2)
makepretty;
hold on;

iAnimal = 5
plot(nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), 'Color', rgb('Blue'))
makepretty;
hold on;
errorbar([1, 2], nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), ...
    nanstd([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]./sqrt(2)), ...
    'Color', rgb('Blue'), 'LineWidth', 2)
makepretty;
hold on;

hold on;
plot(nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), 'Color', rgb('Red'))
makepretty;
hold on;
errorbar([1, 2], nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]), ...
    nanstd([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end -1:end), [1, 3]))]./sqrt(2)), ...
    'Color', rgb('Red'), 'LineWidth', 2)
makepretty;
hold on;

ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
xticks([1, 2, 3])
xticklabels({'Go1', 'NoGo'})
ylim([0, 1])
makepretty;
grid on;

%% phases
figure();
min = [-0.1, 0, 0.1];
for iAnimal = 5:7
    phases = [ones(length(animalsPhase1{iAnimal-4}), 1); ones(length(animalsPhase2{iAnimal-4}), 1) .* 2; ones(length(animalsPhase3{iAnimal-4}), 1) .* 3; ...
        ones(length(animalsPhase4{iAnimal-4}), 1) .* 4; ones(length(animalsPhase5{iAnimal-4}), 1) .* 5];
    p(iAnimal) = plot(phases+min(iAnimal-4), 'Color', animalCol{iAnimal-4});
    makepretty;
    hold on;
end
xlabel('Training day #')
ylabel('Training phase')
legend([p(5), p(6), p(7)], {'JF042', 'JF043', 'JF044'})
makeprettyLarge;

%% n trials
figure();
min = [-0.1, 0, 0.1];
for iAnimal = 5:7
    phases = [sum(an(iAnimal).n_trials(animalsPhase1{iAnimal-4}, :), 2); sum(an(iAnimal).n_trials(animalsPhase2{iAnimal-4}, :), 2); ...
        sum(an(iAnimal).n_trials(animalsPhase3{iAnimal-4}, :), 2); sum(an(iAnimal).n_trials(animalsPhase4{iAnimal-4}, :), 2); sum(an(iAnimal).n_trials(animalsPhase5{iAnimal-4}, :), 2)];
    p(iAnimal) = plot(phases+min(iAnimal-4), 'Color', animalCol{iAnimal-4});
    makepretty;
    hold on;
end
xlabel('Training day')
ylabel('# trials')
legend([p(5), p(6), p(7)], {'JF042', 'JF043', 'JF044'})
makeprettyLarge;
end


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



