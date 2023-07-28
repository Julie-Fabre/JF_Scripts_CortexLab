%clear all;
%clc;
function bhvOut = JF_noGoWorld_behavior_clean(animalsAll)
plotAll = false;
bhvOut = struct;
anBhv = struct;
for iAnimal = 1:size(animalsAll, 2)

    myPaths;
    animal = animalsAll{1, iAnimal}; %this animal


    % get experiments
    protocol = 'stage'; %'location'; %protocol name contains this name
    protocol2 = 'location';
    protocol3 = 'Hack';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments1 = AP_find_experimentsJF(animal, protocol, flexible_name);
    experiments2 = AP_find_experimentsJF(animal, protocol2, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol3, flexible_name);
    experiments(1:size(experiments1, 1)) = experiments1;
    experiments(size(experiments1, 1)+1:size(experiments1, 1)+size(experiments2, 1)) = experiments2;
    experiments(size(experiments1, 1)+size(experiments2, 1)+1:size(experiments1, 1)+size(experiments2, 1)+size(experiments3, 1)) = experiments3;

    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];


    for curr_day = 1:length(experiments)


        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        trialNum = [];

        % If multiple experiments, only use the largest one (usually multiple
        % happens if mess ups first/last one is good)
        for curr_experiment = 1:length(experiment_num)
            try
                experiment = experiment_num(curr_experiment);

                [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
                if ~isempty(block_filename)
                    load(block_filename)

                    trialNum(curr_experiment) = length(block.events.newTrialTimes);
                else
                    trialNum(curr_experiment) = 0;
                end
            catch
            end
        end

        %
        for curr_experiment = find(trialNum == max(trialNum))

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

            % trial side
            if isfield(block.events, 'trialSideValues')
                correctTrialsRight = block.events.hitValues(response_trials) == 1;
                thisSide = 1;
            else
                % get trail side - hardcoded
                if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
                    thisSide = -1;
                else
                    thisSide = 1;
                end

                correctTrialsRight = block.events.feedbackValues(response_trials) == thisSide; %go1 and go2
            end

            % n trials and water amounts
            sum(correctTrialsRight);
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);
            water_amount = block.outputs.rewardValues(1);


            response_trials = 1:length(block.events.endTrialValues);
            clearvars correctTrials


            if isfield(block.events, 'noGoQuiescenceTimes')
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

            % seperate trials by type
            if contains(block.expDef, 'Hack') || contains(block.expDef, 'noGo_stage4') || contains(block.expDef, 'noGo_stage5') || contains(block.expDef, 'path')
                go1Trials = block.events.stimulusTypeValues(response_trials) == 1;
                go2Trials = block.events.stimulusTypeValues(response_trials) == 2;
                noGoTrials = block.events.stimulusTypeValues(response_trials) == 3;
            else
                go1Trials = ones(size(block.events.trialSideValues(response_trials), 2), 1)';
                go2Trials = zeros(size(block.events.trialSideValues(response_trials), 2), 1)';
                noGoTrials = zeros(size(block.events.trialSideValues(response_trials), 2), 1)';
            end

            % get repeat on misses
            if isfield(block.events, 'repeatTrialValues')
                repeatOnMisses = block.events.repeatTrialValues(response_trials) > 0;
            else
                repeatOnMisses = block.events.repeatNumValues(response_trials) > 1;
            end

            % number of correct, non-repeat-on-miss trials
            correctTrials(1) = sum(block.events.hitValues(response_trials) == 1 & go1Trials); % go 1
            correctTrials(2) = sum(block.events.hitValues(response_trials) == 1 & go2Trials); % go 2
            correctTrials(3) = sum(block.events.hitValues(response_trials) == 1 & noGoTrials); % no go

            nTrials(1) = sum(~repeatOnMisses(response_trials) & go1Trials);
            nTrials(2) = sum(~repeatOnMisses(response_trials) & go2Trials);
            nTrials(3) = sum(~repeatOnMisses(response_trials) & noGoTrials);

            stimConditions = [1, 2, 3];

            % number of "go left" non-repeat-on-miss trials
            goLeft(1) = sum(~repeatOnMisses(response_trials) & block.events.hitValues(response_trials) == thisSide & go1Trials);
            goLeft(2) = sum(~repeatOnMisses(response_trials) & block.events.hitValues(response_trials) == thisSide & go2Trials);
            goLeft(3) = sum(~repeatOnMisses(response_trials) & block.events.hitValues(response_trials) == thisSide & noGoTrials);

            % number of "no go" non-repeat-on-miss trials
            noGo(1) = sum(~repeatOnMisses(response_trials) & (block.events.hitValues(response_trials) == 0) & go1Trials);
            noGo(2) = sum(~repeatOnMisses(response_trials) & (block.events.hitValues(response_trials) == 0) & go2Trials);
            noGo(3) = sum(~repeatOnMisses(response_trials) & (block.events.hitValues(response_trials) == 0) & noGoTrials);

            % get wheel aligned to stimulus onset
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1 / wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, wheel_t_resample);

            wheel_smooth_t = 0.05; % seconds
            wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
            wheel_velocity = interp1(conv(wheel_t_resample, [1, 1]/2, 'valid'), ...
                diff(smooth(wheel_values_resample, wheel_smooth_samples)), wheel_t_resample)';
            wheel_thresh = 0.025;


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

            % get wheel move start times
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
            else
                bhv.stim_aligned_wheelGo1{curr_day} = stim_aligned_wheel;
            end
            if sum(go2Trials) > 0
                bhv.stim_aligned_wheelGo2{curr_day} = stim_aligned_wheel(go2Trials, :);

            else
                bhv.stim_aligned_wheelGo2{curr_day} = [];
            end
            if sum(noGoTrials) > 0
                bhv.stim_aligned_wheelNoGo{curr_day} = stim_aligned_wheel(noGoTrials, :);

            else
                bhv.stim_aligned_wheelNoGo{curr_day} = [];
            end
            bhv.correctTrials(curr_day, :) = correctTrials;
            bhv.nTrials(curr_day, :) = nTrials;
            bhv.n_trials(curr_day, :) = n_trials;
            bhv.total_water(curr_day, :) = total_water;
            bhv.water_amount(curr_day, :) = water_amount;
            bhv.conditions(curr_day, :) = stimConditions;
            bhv.goLeft(curr_day, :) = goLeft;
            bhv.noGo(curr_day, :) = noGo;
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

            load_parts.ephys = false;
            load_parts.cam = false;
            isSpikeGlx = false;
            JF_load_experiment;

            expDef_fn = [fileparts(block_filename), filesep, day, '_', ...
                num2str(experiment), '_', animal, '_expDef.m'];
            if ~exist(expDef_fn, 'file')
                error('%s %s: no expDef.m', animal, day)
            end
            expDef_text = fileread(expDef_fn);
            [~, quiescThreshold_txt] = regexp(expDef_text, ...
                'quiescThreshold = (\d*)', 'match', 'tokens');
            quiescThreshold = str2num(quiescThreshold_txt{1}{1});

            % Resetting quiescence period for each trial %qq
            if isfield(block.events, 'trialQuiescenceValues')
                try
                    quiescence_t = block.events(1:n_trials).trialQuiescenceValues;
                catch
                    quiescence_t = block.events.trialQuiescenceValues(1:n_trials);
                end
                iti_t = block.events.trialITIValues;
                t = Timeline.rawDAQTimestamps';
                try
                     quiescence_reset_t = JF_extrap_stimWheel_quiescence(n_trials(end), signals_events, t, block2timeline, timeline2block, block, quiescThreshold);
               
                    [alt_stim_to_move, alt_stimOn_times] = JF_randsampleITI_stimToMove(n_trials(end), signals_events.responseTimes, ...
                    block.paramsValues(1).itiMin, block.paramsValues(1).itiMax, block.paramsValues(1).quiescMin, ...
                    block.paramsValues(1).quiescMax, quiescence_reset_t, t, wheel_starts, wheel_move_stim_idx, quiescence_t, iti_t, signals_events, block);
                    alt_stim_to_move_all = cell2mat(alt_stim_to_move);
                
                catch
                   
                end
                if contains(block.expDef, 'NoGo') || contains(block.expDef, 'stage4') || contains(block.expDef, 'stage5') 

                    alt_stim_to_move_resampled_go1 = datasample(alt_stim_to_move_all(go1Trials(1:end)), 10000);
                    bhv.alt_stim_to_move_resampled(curr_day, 1, :) = alt_stim_to_move_resampled_go1;
                    alt_stim_to_move_resampled_go2 = datasample(alt_stim_to_move_all(go2Trials(1:end)), 10000);
                    bhv.alt_stim_to_move_resampled(curr_day, 2, :) = alt_stim_to_move_resampled_go2;
                    alt_stim_to_move_resampled_noGo = datasample(alt_stim_to_move_all(noGoTrials(1:end)), 10000);
                    bhv.alt_stim_to_move_resampled(curr_day, 3, :) = alt_stim_to_move_resampled_noGo;
                else
                    if ~isempty(alt_stim_to_move_all)
                        alt_stim_to_move_resampled = datasample(alt_stim_to_move_all, 10000);
                        bhv.alt_stim_to_move_resampled(curr_day, 1, :) = alt_stim_to_move_resampled;
                        bhv.alt_stim_to_move_resampled(curr_day, 2, :) = nan(size(alt_stim_to_move_resampled, 1), 1);
                        bhv.alt_stim_to_move_resampled(curr_day, 3, :) = nan(size(alt_stim_to_move_resampled, 1), 1);
                    else
                        bhv.alt_stim_to_move_resampled(curr_day, 1, :) = nan(10000, 1);
                        bhv.alt_stim_to_move_resampled(curr_day, 2, :) = nan(10000, 1);
                        bhv.alt_stim_to_move_resampled(curr_day, 3, :) = nan(10000, 1);
                    end

                end
            else
                bhv.alt_stim_to_move_resampled(curr_day, 1, :) = nan(10000, 1);
                bhv.alt_stim_to_move_resampled(curr_day, 2, :) = nan(10000, 1);
                bhv.alt_stim_to_move_resampled(curr_day, 3, :) = nan(10000, 1);
                %
            end


            keep_day = [keep_day, curr_day];
            expDefNameStart = strfind(block.expDef, 'stage');
            expDefName = block.expDef(expDefNameStart:expDefNameStart+5);


            if contains(block.expDef, 'NoGo') || contains(block.expDef, 'stage4') || contains(block.expDef, 'stage5') || contains(block.expDef, 'path')
                noGoDay = [noGoDay, curr_day];
                if sum(go2Trials) >= 1

                    bhv.movingFracGo2(curr_day, :) = movingFracGo2;
                end
                bhv.stim_to_moveMean(curr_day, 1) = nanmean(stim_to_move(go1Trials(1:end)));
                if any(go2Trials)
                    bhv.stim_to_move{curr_day, 2} = stim_to_move(go2Trials(1:end));
                    bhv.stim_to_move{curr_day, 3} = stim_to_move(noGoTrials(1:end));
                elseif any(noGoTrials)
                    bhv.stim_to_move{curr_day, 2} = nan(size(movingFrac, 2), 1)';
                    bhv.stim_to_move{curr_day, 3} = stim_to_move(noGoTrials(1:end));
                else
                    bhv.stim_to_move{curr_day, 2} = nan(size(movingFrac, 2), 1)';
                    bhv.stim_to_move{curr_day, 3} = nan(size(movingFrac, 2), 1)';
                end
                try
                bhv.stim_to_moveMean(curr_day, 2) = nanmean(stim_to_move(go2Trials(1:end)));
                catch
                end
                bhv.stim_to_moveMean(curr_day, 3) = nanmean(stim_to_move(noGoTrials(1:end)));

                bhv.stim_to_move{curr_day, 1} = stim_to_move(go1Trials(1:end));
                % bhv.stim_to_move{curr_day,2} = stim_to_move(go2Trials(1:end) );
                % bhv.stim_to_move{curr_day,3} = stim_to_move(noGoTrials(1:end) );


                bhv.movingFracGo1(curr_day, :) = movingFracGo1;
                bhv.movingFracNoGo(curr_day, :) = movingFracNoGo;


            else
                %bhv.movingFrac(curr_day, 1) = movingFrac;
                bhv.stim_to_move{curr_day, 1} = stim_to_move;
                bhv.stim_to_move{curr_day, 2} = nan(size(movingFrac, 2), 1)';
                bhv.stim_to_move{curr_day, 3} = nan(size(movingFrac, 2), 1)';


            end
            bhv.stim_rxn_time(curr_day) = stim_rxn_time;
            bhv.stim_rxn_timeSEM(curr_day) = stim_rxn_timeSEM;

            %stims

            keep_day = [keep_day, curr_day];
            expDefNameStart = strfind(block.expDef, 'stage');
            expDefName = block.expDef(expDefNameStart:expDefNameStart+5);
            bhv.expDefName{curr_day} = expDefName;
        end

    end

    keep_day = unique(keep_day);
    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels_temp = cellfun(@(day, protocol) [day(6:end)], ...
        {experiments(keep_day).day}, bhv.expDefName(keep_day), 'uni', false);
    bhvOut(iAnimal). dates = cellfun(@(day, protocol) [day], ...
        {experiments(keep_day).day}, bhv.expDefName(keep_day), 'uni', false);
    %  day_labels = cellfun(@(day, protocol) [protocol, ' ', day(6:end)], ...
    %    {experiments(keep_day).day},bhv.expDefName(keep_day), 'uni', false);
    %
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

    % if ismember(iAnimal, 1:3)
    %     bhv.goLeft(keep_day, 3) = bhv.goLeft(keep_day, 3) -0.2;
    %     bhv.noGo(keep_day, 3) = bhv.noGo(keep_day, 3) +0.2;
    % end
    goL=bhv.goLeft(keep_day, :)./bhv.nTrials(keep_day, :);
    goL(21:end,:) = 1 - goL(21:end,:);
    subplot(232)
    con = bhv.conditions(keep_day(1), :);
    ccc = num2str(con(1, :));
    cc = textscan(ccc, '%s', 'Delimiter', ' ');
    cc = cc{1};
    cc(strcmp('', cc)) = [];
    im = imagesc(1:length(con), 1:size(bhv.correctTrials(keep_day, :), 1), goL);
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
    caxis([0.7, 1])
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

    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.noGo(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    %axis square;
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
    %     for iDay = 1:size(keep_day, 2)
    %         txt = num2str(bhv.trialTime(keep_day(iDay), :)');
    %         text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    %     end
    caxis([0, 1.8])
    makepretty;
    try
        day_idx_go1 = [];
        subplot(236) % last day histogram of reaction times per condition
        clearvars h1
        for iDay = 1:size(bhv.stim_to_move, 1)
            [h1(iDay, 1, :), bins1(iDay, 1, :)] = hist(bhv.stim_to_move{iDay, 1}, 0:0.05:1.8);
            %day_idx_go1=[day_idx_go1; ones(size(bhv.stim_to_move{iDay,1},2),1)*iDay];
            n_rxn_altsample = 1000;
            rxn_measured_med(iDay, 1) = nanmedian(bhv.stim_to_move{iDay, 1});
            rxn_alt_med(iDay, 1) = nanmedian(datasample(bhv.stim_to_move{iDay, 1}, n_rxn_altsample));
            [h1(iDay, 2, :), bins1(iDay, 2, :)] = hist(bhv.stim_to_move{iDay, 2}, 0:0.05:1.8);
            if ~any(h1(iDay, 2, :))
                h1(iDay, 2, :) = nan(1, 37);
            end
            [h1(iDay, 3, :), bins1(iDay, 3, :)] = hist(bhv.stim_to_move{iDay, 3}, 0:0.05:1.8);
            if ~any(h1(iDay, 3, :))
                h1(iDay, 3, :) = nan(1, 37);
            end
        end
        if any(h1(end, 3, :))

            im = imagesc([squeeze(h1(:, 1, :)), squeeze(h1(:, 2, :)), squeeze(h1(:, 3, :))]);
            set(im, 'AlphaData', ~isnan(get(im, 'CData')));
            set(gca, 'color', [0.5, 0.5, 0.5]);
            xticks([1:18.5:37 * 3 + 1])
            xticklabels(strtrim(cellstr(num2str([0:0.9:1.8 * 3]'))'))

        else
            im = imagesc(squeeze(bins1(iDay, 1, :))', [], squeeze(h1(:, 1, :)));
        end

        colormap(brewermap([], '*RdBu'));
        caxis([- max(abs(max(max(im.CData)))), max(abs(max(max(im.CData))))])
    catch
    end
    ylabel('training day')
    xlabel('reaction time (s)')

    subplot(235)
    stim_to_move_day_median = cellfun(@(x) median(x, 'all', 'omitnan'), bhv.stim_to_move(:, 1));
    % stim_to_move_day_std = cellfun(@(x) median(x-median(x, 'all', 'omitnan'), 'all', 'omitnan'), bhv.stim_to_move(:, 1));
    semilogy(1:length(bhv.stim_to_move), stim_to_move_day_median, 'Color', 'b');
    hold on;
    x = 1:length(bhv.stim_to_move);
%    lo = stim_to_move_day_median - stim_to_move_day_std;
 %   hi = stim_to_move_day_median + stim_to_move_day_std;
%    lo(lo < 0) = 0.0001;
 %   hi(hi < 0) = 0.0001;
%    patch([x, x(end:-1:1), x(1)], [lo; hi(end:-1:1); lo(1)], rgb('Blue'), 'FaceAlpha', .3);


    stim_to_move_alt_day_median = nanmedian(bhv.alt_stim_to_move_resampled(:, 1, :), 3);
    stim_to_move_alt_day_ci = squeeze(prctile(bhv.alt_stim_to_move_resampled(:, 1, :), [5, 95], 3));

    semilogy(1:length(bhv.stim_to_move), stim_to_move_alt_day_median, 'Color', rgb('Grey'));
    x = 1:length(bhv.stim_to_move);
    lo = stim_to_move_alt_day_ci(:, 1);
    hi = stim_to_move_alt_day_ci(:, 2);
    lo(lo < 0, :) = 0.0001;
    hi(hi < 0, :) = 0.0001;
    x(isnan(lo(:, 1))) = [];
    hi(isnan(lo(:, 1)), :) = [];
    lo(isnan(lo(:, 1)), :) = [];
    %patch([x, x(end:-1:1), x(1)], [lo; hi(end:-1:1); lo(1)], rgb('Gray'),'FaceAlpha',.3);
    % shading shows 95% confidence intervals from the null distribution
    % median absolute deviation across mice
    makepretty;
    xlabel('training day')
    ylabel('reaction time (s)')
    ylim([0.1, 4])


    bhvOut(iAnimal).protocol_idx = protocol_idx;
    bhvOut(iAnimal).day_labels = day_labels;
    bhvOut(iAnimal).binBorders = binBorders;
    bhvOut(iAnimal).goLeft = bhv.goLeft(:, :);
    bhvOut(iAnimal).nTrials = bhv.nTrials(:, :);
    bhvOut(iAnimal).noGo = bhv.noGo(:, :);
    bhvOut(iAnimal).stim_to_move = stim_to_move;
    bhvOut(iAnimal).noGoDay = noGoDay;
    bhvOut(iAnimal).stim_aligned_wheelGo1 = bhv.stim_aligned_wheelGo1;
    bhvOut(iAnimal).stim_aligned_wheelGo2 = bhv.stim_aligned_wheelGo2;
    bhvOut(iAnimal).stim_aligned_wheelNoGo = bhv.stim_aligned_wheelNoGo;
    bhvOut(iAnimal).repeatOnMisses = bhv.repeatOnMisses;
    bhvOut(iAnimal).response_trials = bhv.response_trials;
    bhvOut(iAnimal).hitValues = bhv.hitValues;
    bhvOut(iAnimal).go1trials = bhv.go1trials;
    bhvOut(iAnimal).noGotrials = bhv.noGotrials;
    bhvOut(iAnimal).stim_to_move = bhv.stim_to_move;
    anBhv(iAnimal).noGo = bhv.noGo;
    anBhv(iAnimal).nTrials = bhv.nTrials;
    anBhv(iAnimal).goLeft = bhv.goLeft;
    anBhv(iAnimal).noGoDay = noGoDay;
    anBhv(iAnimal).stim_rxn_time = bhv.stim_rxn_time;
    anBhv(iAnimal).n_trials_all = bhv.n_trials;
    clearvars bhv

end

%% summary: % no go per stim per mouse on last day only


for iAnimal = 1:length(animalsAll)
    try
        %subplot(3, 1, iAnimal)
        figure();
        transparencyValues = 0:1 / length(anBhv(iAnimal).noGoDay):1;
        for iDay = 1:size(anBhv(iAnimal).noGoDay, 2)
            if iDay == size(anBhv(iAnimal).noGoDay, 2)
                col1 = rgb('DarkBlue');
                col2 = rgb('DarkRed');

            else
                col1 = 'b';
                col2 = 'r';
            end
            scatter([1, 2, 3], anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2])./anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2]), [], col1, 'filled', ...
                'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
            p = plot([1, 2, 3], anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2])./anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2]), 'Color', col1);
            p.Color(4) = transparencyValues(iDay+1);
            hold on;
            makepretty;
            scatter([1, 2, 3], anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2])./anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2]), [], col2, 'filled', ...
                'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
            p = plot([1, 2, 3], anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2])./anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(iDay), [1, 3, 2]), 'Color', col2);
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
        anBhv(iAnimal).noGoDay = unique(anBhv(iAnimal).noGoDay);
        subplot(round(length(animalsAll)/2), round(length(animalsAll)/2), iAnimal)

        transparencyValues = 0:1 / length(anBhv(iAnimal).noGoDay):1;
        %for iDay = 1:size(an(iAnimal).noGoDay, 2)
        col1 = rgb('DarkBlue');
        col2 = rgb('DarkRed');


        scatter([1, 2, 3], nanmean(anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), [], col1, 'filled');
        hold on;
        p = plot([1, 2, 3], nanmean(anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color', col1);

        hold on;
        %errorbar([1, 2, 3],nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])),...
        %    nanstd(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])./an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color',col1)
        makepretty;
        scatter([1, 2, 3], nanmean(anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), [], col2, 'filled')
        p = plot([1, 2, 3], nanmean(anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), 'Color', col2);

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

clearvars meanGoLeftAcrossAnimals
clearvars meanNoGoAcrossAnimals
figure();
for iAnimal = 1:length(animalsAll)
    try
        transparencyValues = 0:1 / length(anBhv(iAnimal).noGoDay):1;
        %for iDay = 1:size(an(iAnimal).noGoDay, 2)
        col1 = rgb('DarkBlue');
        col2 = rgb('DarkRed');
        if ismember(iAnimal, 1:4)
            meanGoLeftAcrossAnimals(iAnimal, :) = nanmean(anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])) ./ nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
            meanNoGoAcrossAnimals(iAnimal, :) = nanmean(anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])) ./ nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
        else
            meanGoLeftAcrossAnimals(iAnimal, :) = 1 - nanmean(anBhv(iAnimal).goLeft(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])) ./ nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
            meanNoGoAcrossAnimals(iAnimal, :) = nanmean(anBhv(iAnimal).noGo(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])) ./ nanmean(anBhv(iAnimal).nTrials(anBhv(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));

        end
    catch
    end
end
try
    scatter([1, 2, 3], nanmean(meanGoLeftAcrossAnimals), [], col1, 'filled');
    hold on;
    p = plot([1, 2, 3], nanmean(meanGoLeftAcrossAnimals), 'Color', col1);
    hold on;
    errorbar([1, 2, 3], nanmean(meanGoLeftAcrossAnimals), ...
        nanstd(meanGoLeftAcrossAnimals)./length(animalsAll), 'Color', col1)
    makepretty;

    scatter([1, 2, 3], nanmean(meanNoGoAcrossAnimals), [], col2, 'filled')
    p = plot([1, 2, 3], nanmean(meanNoGoAcrossAnimals), 'Color', col2);

    errorbar([1, 2, 3], nanmean(meanNoGoAcrossAnimals), ...
        nanstd(meanNoGoAcrossAnimals)./length(animalsAll), 'Color', col2)
    makepretty;
    axis square

    xlabel('image type')
    ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
    xticks([1, 2, 3])
    xticklabels({'Go1', 'NoGo', 'Go2'})
    ylim([0, 1])
    dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    legend(dummyh, ['mean +/- s.e. across ', num2str(length(animalsAll)), ' mice'], 'Box', 'off')
    makepretty;
    grid on;
catch
end

end