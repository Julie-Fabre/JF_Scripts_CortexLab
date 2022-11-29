%clear all;
%clc;
function bhvOut = JF_noGoWorld_behavior(animalsAll)

bhvOut = struct;
for iAnimal = 1:size(animalsAll, 2)


    animal = animalsAll{1, iAnimal}; %this animal
    if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
        thisSide = -1;
    else
        thisSide = 1;
    end
    protocol = 'stage'; %'location'; %protocol name contains this name
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, flexible_name);
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];

    for curr_day = 1:length(experiments)


        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        trNum = [];
        % If multiple experiments, only use the largest one (usually multiple
        % happens if mess ups first/last one is good)
        for curr_experiment = 1:length(experiment_num)
            try
                experiment = experiment_num(curr_experiment);

                [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
                if ~isempty(block_filename)
                    load(block_filename)

                    trNum(curr_experiment) = length(block.events.newTrialTimes);
                else
                    trNum(curr_experiment) = 0;
                end
            catch
            end

            %
            for curr_experiment = find(trNum == max(trNum))
                try
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
                        correctTrialsRight = block.events.hitValues(response_trials) == 1;
                    else
                        correctTrialsRight = block.events.feedbackValues(response_trials) == thisSide; %go1 and go2
                    end
                    sum(correctTrialsRight);
                    n_trials = length(block.paramsValues);
                    total_water = sum(block.outputs.rewardValues);
                    water_amount = block.outputs.rewardValues(1);


                    response_trials = 1:length(block.events.endTrialValues);
                    clearvars correctTrials
                    iImg = 1;


                    if isfield(block.events, 'trialSideValues')
                        %leftGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == -1));
                        rightGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == thisSide));
                        %leftT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == -1));
                        rightT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == thisSide));
                    elseif isfield(block.events, 'noGoQuiescenceTimes')
                        block.events.trialSideValues(response_trials) = thisSide;
                        hitTimesPossib = [block.events.feedbackTimes, block.events.noGoQuiescenceTimes];
                        hitValuesUnsorted = [block.events.feedbackValues, 0 * ones(size(block.events.noGoQuiescenceTimes, 2), thisSide)'];
                        [sortedV, sortedIdx] = sort(hitTimesPossib);
                        block.events.hitValues(response_trials) = hitValuesUnsorted(sortedIdx(response_trials));
                    else
                        block.events.trialSideValues(response_trials) = thisSide;
                        ff = block.events.responseValues;

                        block.events.hitValues(response_trials) = ff(response_trials);
                    end
                    if contains(block.expDef, 'goNoGo')
                        go1Trials = block.events.trialTypeValues(response_trials) == 1;
                        go2Trials = block.events.trialTypeValues(response_trials) == 2;
                        noGoTrials = block.events.trialTypeValues(response_trials) == 3;
                    elseif contains(block.expDef, 'Hack') || contains(block.expDef, 'noGo_stage4') || contains(block.expDef, 'noGo_stage5')
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
                    correctTrials((iImg - 1)*3+1) = sum( ...
                        block.events.hitValues(response_trials) == 1 & go1Trials ...
                        ); % go 1
                    correctTrials((iImg - 1)*3+2) = sum( ...
                        block.events.hitValues(response_trials) == 1 & go2Trials ...
                        ); % go 2
                    correctTrials((iImg - 1)*3+3) = sum( ...
                        block.events.hitValues(response_trials) == 1 & noGoTrials ...
                        ); % no go

                    nTrials((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                        go1Trials ...
                        );
                    nTrials((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                        go2Trials ...
                        );
                    nTrials((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                        noGoTrials ...
                        );

                    conditi = [1, 2, 3];


                    if contains(block.expDef, 'Hack') || contains(block.expDef, 'noGo_stage4') || contains(block.expDef, 'noGo_stage5')
                        if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
                            val = 1;
                        else
                            val = -1;
                        end

                    else
                        if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
                            val = -1;
                        else

                            val = 1;
                        end
                    end
                    goLeft((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go1Trials ...
                        );
                    goLeft((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go2Trials ...
                        );
                    goLeft((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & noGoTrials ...
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
                        ...
                        go1Trials ...
                        ));
                    trialTime((iImg - 1)*3+2) = nanmean(timeThis(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                        ...
                        go2Trials ...
                        ));
                    trialTime((iImg - 1)*3+3) = nanmean(timeThis(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                        ...
                        noGoTrials ...
                        ));


                    noGo((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        (block.events.hitValues(response_trials) == 0) & go1Trials ...
                        );
                    noGo((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        (block.events.hitValues(response_trials) == 0) & go2Trials ...
                        );
                    noGo((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        (block.events.hitValues(response_trials) == 0) & noGoTrials ...
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
                    surround_time = [-0.5, 2];
                    surround_sample_rate = 1 / 0.001; % (match this to framerate)
                    surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
                    if isfield(block.events, 'stimOnTimes')
                        pull_times = bsxfun(@plus, block.events.stimOnTimes', surround_time_points);
                    else
                        pull_times = bsxfun(@plus, block.events.stimulusOnTimes', surround_time_points);
                    end

                    stim_aligned_wheel = interp1(wheel_t_resample, ...
                        wheel_velocity, pull_times);


                    % (set a threshold in speed and time for wheel movement)
                    thresh_displacement = 0.025; %% QQ SHOULD BE 0.025
                    time_over_thresh = 0.05; % ms over velocity threshold to count
                    samples_over_thresh = time_over_thresh .* surround_sample_rate;
                    wheel_over_thresh_fullconv = convn( ...
                        abs(stim_aligned_wheel) > thresh_displacement, ...
                        ones(1, samples_over_thresh)) >= samples_over_thresh;
                    wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

                    [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
                    wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
                    wheel_move_time(~move_trial) = NaN;

                    if isfield(block.events, 'stimOnTimes')
                        stim_to_move = padarray(wheel_move_time'-block.events.stimOnTimes, ...
                            [n_trials - length(block.events.stimOnTimes), 0], NaN, 'post');

                        trial_wheel_starts = arrayfun(@(x) ...
                            wheel_starts(find(wheel_starts > block.events.stimOnTimes(x), 1)), ...
                            response_trials);

                        trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
                        stim_rxn_time = nanmedian(trial_move_t);
                        stim_rxn_timeSEM = nanstd(trial_move_t) ./ sqrt(numel(trial_move_t));
                    else
                        try
                            stim_to_move = padarray(wheel_move_time'-block.events.stimulusOnTimes, ...
                                [n_trials - length(block.events.stimulusOnTimes), 0], NaN, 'post');
                            trial_wheel_starts = arrayfun(@(x) ...
                                wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                                response_trials(1:end));
                            trial_wheel_move_start = arrayfun(@(x) ...
                                find(wheel_t_resample > block.events.stimulusOnTimes(x), 1), ...
                                response_trials(1:end));
                            trial_wheel_move = nan(size(response_trials(1:end), 2), 1600); %about 1.6 seconds
                            for iTrial = 1:size(response_trials(1:end), 2)
                                trial_wheel_move(iTrial, :) = wheel_velocity(trial_wheel_move_start(iTrial): ...
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
                            trial_wheel_move = nan(size(response_trials(1:end), 2), 1600);
                        end
                    end


                    movingRN = abs(wheel_velocity) > 0.022;

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
                    if sum(noGoTrials) > 0
                        bhv.stim_aligned_wheelGo1{curr_day} = stim_aligned_wheel(go1Trials, :);
                        bhv.stim_aligned_wheelGo1_inco{curr_day} = stim_aligned_wheel(go1Trials & block.events.hitValues(response_trials) ~= val, :);
                        bhv.stim_aligned_wheelGo1_corr{curr_day} = stim_aligned_wheel(go1Trials & block.events.hitValues(response_trials) == val, :);
                    else
                        bhv.stim_aligned_wheelGo1{curr_day} = stim_aligned_wheel;
                        bhv.stim_aligned_wheelGo1_corr{curr_day} = [];
                        bhv.stim_aligned_wheelGo1_inco{curr_day} = [];
                    end
                    if sum(go2Trials) > 0
                        bhv.stim_aligned_wheelGo2{curr_day} = stim_aligned_wheel(go2Trials, :);
                        bhv.stim_aligned_wheelGo2_inco{curr_day} = stim_aligned_wheel(go2Trials & block.events.hitValues(response_trials) ~= val, :);
                        bhv.stim_aligned_wheelGo2_corr{curr_day} = stim_aligned_wheel(go2Trials & block.events.hitValues(response_trials) == val, :);

                    else
                        bhv.stim_aligned_wheelGo2{curr_day} = [];
                        bhv.stim_aligned_wheelGo2_corr{curr_day} = [];
                        bhv.stim_aligned_wheelGo2_inco{curr_day} = [];
                    end
                    if sum(noGoTrials) > 0
                        bhv.stim_aligned_wheelNoGo{curr_day} = stim_aligned_wheel(noGoTrials, :);
                        bhv.stim_aligned_wheelNoGo_corr{curr_day} = stim_aligned_wheel(noGoTrials & block.events.hitValues(response_trials) == 0, :);
                        bhv.stim_aligned_wheelNoGo_inco{curr_day} = stim_aligned_wheel(noGoTrials & block.events.hitValues(response_trials) ~= 0, :);

                    else
                        bhv.stim_aligned_wheelNoGo{curr_day} = [];
                        bhv.stim_aligned_wheelNoGo_corr{curr_day} = [];
                        bhv.stim_aligned_wheelNoGo_inco{curr_day} = [];
                    end
                    bhv.correctTrials(curr_day, :) = correctTrials;
                    bhv.nTrials(curr_day, :) = nTrials;
                    bhv.n_trials(curr_day, :) = n_trials;
                    bhv.total_water(curr_day, :) = total_water;
                    bhv.water_amount(curr_day, :) = water_amount;
                    bhv.conditions(curr_day, :) = conditi;
                    bhv.goLeft(curr_day, :) = goLeft;
                    bhv.noGo(curr_day, :) = noGo;
                    bhv.goLeft1(curr_day, :) = goLeft1;
                    bhv.noGo1(curr_day, :) = noGo1;
                    bhv.goLeft2(curr_day, :) = goLeft2;
                    bhv.noGo2(curr_day, :) = noGo2;
                    bhv.nTrials1(curr_day, :) = nTrials1;
                    bhv.nTrials2(curr_day, :) = nTrials2;
                    bhv.trialTime(curr_day, :) = trialTime;
                    if size(stim_to_move, 1) == 2
                        stim_to_move = stim_to_move(1, :);
                    end
                    bhv.stim_to_move{curr_day, 1} = stim_to_move;
                    bhv.stim_to_moveMean(curr_day, 1) = nanmean(stim_to_move);
                    bhv.stim_to_moveMean(curr_day, 2) = NaN;
                    bhv.stim_to_moveMean(curr_day, 3) = NaN;
                    bhv.stim_aligned_wheel{curr_day} = stim_aligned_wheel;
                    bhv.repeatOnMisses{curr_day} = repeatOnMisses;
                    bhv.response_trials{curr_day} = response_trials;
                    bhv.hitValues{curr_day} = block.events.hitValues;
                    bhv.go1trials{curr_day} = go1Trials;
                    bhv.noGotrials{curr_day} = noGoTrials;

                    keep_day = [keep_day, curr_day];
                    expDefNameStart = strfind(block.expDef, 'stage');
                    expDefName = block.expDef(expDefNameStart:expDefNameStart+5);
                    bhv.expDefName{curr_day} = expDefName;

                    if contains(block.expDef, 'NoGo') || contains(block.expDef, 'stage4') || contains(block.expDef, 'stage5')
                        noGoDay = [noGoDay, curr_day];
                        if sum(go2Trials) >= 1

                            bhv.movingFracGo2(curr_day, :) = movingFracGo2;
                        end
                        bhv.stim_to_moveMean(curr_day, 1) = nanmean(stim_to_move(go1Trials(1:end)));
                        bhv.stim_to_moveMean(curr_day, 2) = nanmean(stim_to_move(go2Trials(1:end)));
                        bhv.stim_to_moveMean(curr_day, 3) = nanmean(stim_to_move(noGoTrials(1:end)));

                        bhv.stim_to_move{curr_day, 1} = stim_to_move(go1Trials(1:end));
                        bhv.stim_to_move{curr_day, 2} = stim_to_move(go2Trials(1:end));
                        bhv.stim_to_move{curr_day, 3} = stim_to_move(noGoTrials(1:end));


                        bhv.movingFracGo1(curr_day, :) = movingFracGo1;
                        bhv.movingFracNoGo(curr_day, :) = movingFracNoGo;
                        bhv.trialMoveNoGoCorrect{curr_day, :, :} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                        bhv.trialMoveNoGoIncorrect{curr_day, :, :} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                        bhv.trialMoveGo1Correct{curr_day, :, :} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                        bhv.trialMoveGo2Correct{curr_day, :, :} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                        bhv.trialMoveGo1Incorrect{curr_day, :, :} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                        bhv.trialMoveGo2Incorrect{curr_day, :, :} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                    else
                        bhv.movingFrac(curr_day, :) = movingFrac;


                    end
                    bhv.stim_rxn_time(curr_day) = stim_rxn_time;
                    bhv.stim_rxn_timeSEM(curr_day) = stim_rxn_timeSEM;

                    %stims
                catch
                end

            end
        end
    end
    keep_day = unique(keep_day);
    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels_temp = cellfun(@(day, protocol) [day(6:end)], ...
        {experiments(keep_day).day}, bhv.expDefName(keep_day), 'uni', false);
    bhvOut(iAnimal). dates = cellfun(@(day, protocol) [day], ...
        {experiments(keep_day).day}, bhv.expDefName(keep_day), 'uni', false);

    tempLabel = "\\color[rgb]{%s}";


    [unique_protocols, ~, protocol_idx] = unique(bhv.expDefName(keep_day));
    protocol_col = hsv(length(unique_protocols));
    for iDay = 1:length(protocol_idx)
        thisColor = sprintf(tempLabel, num2str(protocol_col(protocol_idx(iDay), :)));
        day_labels{iDay} = append(thisColor, day_labels_temp{iDay});
    end
    figure('Name', animal);
    subplot(231)
    yyaxis left

    scatter(day_num, bhv.n_trials(keep_day, :));
    hold on;
    plot(day_num, bhv.n_trials(keep_day, :), 'linewidth', 2);
    plot(day_num, bhv.water_amount(keep_day, :)*100, 'linewidth', 2, 'color', [1.0000, 0.6445, 0]);
    ylabel('Trials');
    yyaxis right

    scatter(day_num, bhv.total_water(keep_day, :), 'linewidth', 2);
    hold on;
    plot(day_num, bhv.total_water(keep_day, :), 'linewidth', 2);

    ax = gca;
    hold on;

    ylabel('Total water (ul)');
    xlabel('Session');
    set(gca, 'XTick', day_num);
    set(gca, 'XTickLabel', day_labels);
    set(gca, 'XTickLabelRotation', 45);
    %set(gca, 'Color',
    makepretty;


    subplot(232)
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
    caxis([0, 1])
    makepretty;

    subplot(233)
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


    caxis([0, 1])
    makepretty;

    subplot(234)
    im = imagesc(1:length(con), 1:size(nanmean(bhv.stim_to_moveMean(keep_day, :)), 1), bhv.stim_to_moveMean(keep_day, :));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'Mean RT');
    xlabel('image type');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);

    caxis([0, 1.8])
    makepretty;
    try
        subplot(236) % last day histogram of reaction times per condition
        [h1, bins1] = hist(bhv.stim_to_move{end, 1}, 0:0.05:1.8);
        stairs(bins1, h1, 'b', 'LineWidth', 2);
        hold on;
        [h2, bins2] = hist(bhv.stim_to_move{end, 2}, 0:0.05:1.8);
        stairs(bins2, h2, 'r', 'LineWidth', 2);
        [h3, bins3] = hist(bhv.stim_to_move{end, 3}, 0:0.05:1.8);
        stairs(bins3, h3, 'k', 'LineWidth', 2);
        yy = ylim;
        line([0.25, 0.25], [yy(1), yy(2)], 'Color', rgb('Green'));
        ylim([yy(1), yy(2)])
        xlim([0, 1.8])
        legend({['fraction > 0.25 =', num2str(sum(h1(find(bins1 >= 0.25, 1, 'first'):end))./sum(h1))], ...
            ['fraction > 0.25 =', num2str(sum(h2(find(bins2 >= 0.25, 1, 'first'):end))./sum(h2))], ...
            ['fraction > 0.25 =', num2str(sum(h3(find(bins3 >= 0.25, 1, 'first'):end))./sum(h3))]})
        makepretty;
    catch
    end
    subplot(235)

    transparencyValues = 0:1 / length(noGoDay):1;
    bhvOut(iAnimal).protocol_idx = protocol_idx;
    bhvOut(iAnimal).day_labels = day_labels;
    bhvOut(iAnimal).binBorders = binBorders;
    bhvOut(iAnimal).movingFracGo1 = bhv.movingFrac;
    bhvOut(iAnimal).goLeft = bhv.goLeft(:, :);
    bhvOut(iAnimal).nTrials = bhv.nTrials(:, :);
    bhvOut(iAnimal).noGo = bhv.noGo(:, :);
    bhvOut(iAnimal).stim_to_move = stim_to_move;
    bhvOut(iAnimal).noGoDay = noGoDay;
    bhvOut(iAnimal).stim_aligned_wheelGo1 = bhv.stim_aligned_wheelGo1;
    bhvOut(iAnimal).stim_aligned_wheelGo2 = bhv.stim_aligned_wheelGo2;
    bhvOut(iAnimal).stim_aligned_wheelNoGo = bhv.stim_aligned_wheelNoGo;
    bhvOut(iAnimal).stim_aligned_wheelGo1_corr = bhv.stim_aligned_wheelGo1_corr;
    bhvOut(iAnimal).stim_aligned_wheelGo2_corr = bhv.stim_aligned_wheelGo2_corr;
    bhvOut(iAnimal).stim_aligned_wheelNoGo_corr = bhv.stim_aligned_wheelNoGo_corr;
    bhvOut(iAnimal).stim_aligned_wheelGo1_inco = bhv.stim_aligned_wheelGo1_inco;
    bhvOut(iAnimal).stim_aligned_wheelGo2_inco = bhv.stim_aligned_wheelGo2_inco;
    bhvOut(iAnimal).stim_aligned_wheelNoGo_inco = bhv.stim_aligned_wheelNoGo_inco;
    bhvOut(iAnimal).repeatOnMisses = bhv.repeatOnMisses;
    bhvOut(iAnimal).response_trials = bhv.response_trials;
    bhvOut(iAnimal).hitValues = bhv.hitValues;
    bhvOut(iAnimal).go1trials = bhv.go1trials;
    bhvOut(iAnimal).noGotrials = bhv.noGotrials;
    bhvOut(iAnimal).stim_to_move = bhv.stim_to_move;

    if length(noGoDay) >= 1
        bhvOut(iAnimal).movingFracGo1(noGoDay, :) = bhv.movingFracGo1(noGoDay, :);
        bhvOut(iAnimal).movingFracNoGo = bhv.movingFracNoGo;
    else
        bhvOut(iAnimal).movingFracNoGo = zeros(size(bhv.movingFrac, 1), size(bhv.movingFrac, 2));
    end

    if length(noGoDay) == 0


        transparencyValues = 0:1 / length(keep_day):1;
        cla;
        for iDay = length(keep_day)
            p1 = plot(binBorders(1:end-1), bhv.movingFrac(keep_day(iDay), :), 'b');
            p1.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
        end
        ylim([0, 1])
    else
        for iDay = length(noGoDay)
            p1 = plot(binBorders(1:end-1), bhv.movingFracGo1(noGoDay(iDay), :), 'b');
            p1.Color(4) = transparencyValues(iDay+1);
            makepretty;
            hold on;
            ylim([0, 1])

        end
        if isfield(bhv, 'movingFracGo2')
            try
                p2 = plot(binBorders(1:end-1), bhv.movingFracGo2(noGoDay(iDay), :), 'r');
                p2.Color(4) = transparencyValues(iDay+1);
                makepretty;
                hold on;
                legend([p1, p2], {'Go1', 'Go2'})
            catch
            end
        end
        ylim([0, 1])
    end
    xlim([binBorders(1), binBorders(end)])

    xlabel('time from stim onset (s)')
    ylabel('fraction moving')
    makepretty;

    %subplot(336)
    if length(noGoDay) == 0
    else
        transparencyValues = 0:1 / length(noGoDay):1;
        for iDay = length(noGoDay)

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

    for iNoGoday = 1:size(noGoDay, 2)
        try
            meanGo1Correct(iNoGoday, :) = nanmean(abs(bhv.trialMoveGo1Correct{noGoDay(iNoGoday), :, :}));
            meanGo2Correct(iNoGoday, :) = nanmean(abs(bhv.trialMoveGo2Correct{noGoDay(iNoGoday), :, :}));
            meanNoGoCorrect(iNoGoday, :) = nanmean(abs(bhv.trialMoveNoGoCorrect{noGoDay(iNoGoday), :, :}));
            meanGo1Incorrect(iNoGoday, :) = nanmean(abs(bhv.trialMoveGo1Incorrect{noGoDay(iNoGoday), :, :}));
            meanGo2Incorrect(iNoGoday, :) = nanmean(abs(bhv.trialMoveGo2Incorrect{noGoDay(iNoGoday), :, :}));
            meanNoGoIncorrect(iNoGoday, :) = nanmean(abs(bhv.trialMoveNoGoIncorrect{noGoDay(iNoGoday), :, :}));
        catch
            try
                meanGo1Correct(iNoGoday, :) = nan(1600, 1);
            catch
                meanGo1Correct(iNoGoday, :) = NaN;
            end
            meanGo2Correct(iNoGoday, :) = NaN;
            meanNoGoCorrect(iNoGoday, :) = nan(1600, 1);
            meanGo1Incorrect(iNoGoday, :) = nan(1600, 1);
            meanGo2Incorrect(iNoGoday, :) = nan(1600, 1);
            meanNoGoIncorrect(iNoGoday, :) = nan(1600, 1);


        end
    end

    an(iAnimal).noGo = bhv.noGo;
    an(iAnimal).nTrials = bhv.nTrials;
    an(iAnimal).goLeft = bhv.goLeft;
    an(iAnimal).noGoDay = noGoDay;

    an(iAnimal).stim_rxn_time = bhv.stim_rxn_time;
    an(iAnimal).movingFrac = bhv.movingFrac;
    an(iAnimal).n_trials_all = bhv.n_trials;
    clearvars bhv

end

%% summary: % no go per stim per mouse on last day only


for iAnimal = 1:length(animalsAll)
    try
        %subplot(3, 1, iAnimal)
        figure();
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
    catch
    end
end
figure();
clf
clearvars min
try
    nDays = min(1, length(animalsAll));
catch

    nDays = 15;
end

for iAnimal = 1:length(animalsAll)
    try
        an(iAnimal).noGoDay = unique(an(iAnimal).noGoDay);
        subplot(round(length(animalsAll)/2), round(length(animalsAll)/2), iAnimal)

        transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
        %for iDay = 1:size(an(iAnimal).noGoDay, 2)
        col1 = rgb('DarkBlue');
        col2 = rgb('DarkRed');


        scatter([1, 2, 3], nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), [], col1, 'filled');
        hold on;
        p = plot([1, 2, 3], nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color', col1);

        hold on;
        %errorbar([1, 2, 3],nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])),...
        %    nanstd(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color',col1)
        makepretty;
        scatter([1, 2, 3], nanmean(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), [], col2, 'filled')
        p = plot([1, 2, 3], nanmean(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color', col2);

        % errorbar([1, 2, 3],nanmean(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])),...
        %    nanstd(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color',col2)

        makepretty;
        axis square
        %end
        if iAnimal == 1
            xlabel('image type')
            ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
        end
        xticks([1, 2, 3])
        xticklabels({'Go1', 'NoGo', 'Go2'})
        title(animalsAll{iAnimal})
        makepretty;
        ylim([0, 1])

        grid on;
    catch
    end
end


end