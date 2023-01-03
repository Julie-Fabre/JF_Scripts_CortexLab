 expDef = strrep(block.expDef, '\', '/'); % for windows to UNIX paths
    [~, expDef] = fileparts(expDef);
    switch expDef
        case {'AP_stimWheelRight', 'AP_stimWheelLeft', 'AP_stimWheelLeftReverse'}
            % Hit/miss recorded for previous trial, circshift to align
            signals_events.hitValues = circshift(signals_events.hitValues, [0, -1]);
            signals_events.missValues = circshift(signals_events.missValues, [0, -1]);

            % Get number of completed trials (if uncompleted last trial)
            n_trials = length(signals_events.endTrialTimes);

            % Get stim on times by closest photodiode flip
            [~, closest_stimOn_photodiode] = ...
                arrayfun(@(x) min(abs(signals_events.stimOnTimes(x)- ...
                photodiode_flip_times)), ...
                1:n_trials);
            stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);

            % Check that the stim times aren't off by a certain threshold
            % (skip the first one - that's usually delayed a little)
            stim_time_offset_thresh = 0.05;
            if any(abs(stimOn_times(2:end)-signals_events.stimOnTimes(2:n_trials)') >= ...
                    stim_time_offset_thresh)
                figure;
                plot(stimOn_times-signals_events.stimOnTimes(1:n_trials)', '.k')
                line(xlim, repmat(stim_time_offset_thresh, 2, 1), 'color', 'r');
                line(xlim, repmat(-stim_time_offset_thresh, 2, 1), 'color', 'r');
                warning('Stim signals/photodiode offset over threshold');
                xlabel('Stim number');
                ylabel('Photodiode - signals stim time');
                title([animal, ' ', day, ' ', num2str(experiment)]);
            end
            rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
            % (this is a very strange hack to overcome a problem in the rotary
            % encoder that's known in the lab and was put on the wiki)
            wheel_position = Timeline.rawDAQData(:, rotaryEncoder_idx);
            wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;
            [wheel_velocity, wheel_move] = AP_parse_wheel(wheel_position, Timeline.hw.daqSampleRate);


            % Get wheel movement on/offsets
            wheel_starts = Timeline.rawDAQTimestamps(diff([0; wheel_move]) == 1)';
            wheel_stops = Timeline.rawDAQTimestamps(diff([wheel_move; 0]) == -1)';

            % (stim move: first move after stim)
            % (give this a little leeway, sometimes movement starts early but
            % stim comes on anyway)
            stim_leeway = 0.1;
            wheel_move_stim_idx = ...
                arrayfun(@(stim) find(wheel_starts > stim-stim_leeway, 1, 'first'), ...
                stimOn_times);
            wheel_types = ones(size(wheel_starts, 1), 1);
            wheel_types(wheel_move_stim_idx) = 2;


            % (response move: last move start before response signal)
            wheel_move_response_idx = ...
                arrayfun(@(response) find(wheel_starts <= response, 1, 'last'), ...
                signals_events.responseTimes(1:n_trials)');

            % (iti move: move start with no stim on screen)
            stimOff_times = signals_events.stimOffTimes';
            stimOn_epochs = logical(interp1([0; stimOn_times; stimOff_times], ...
                [0; ones(size(stimOn_times)); zeros(size(stimOff_times))], ...
                Timeline.rawDAQTimestamps', 'previous', 'extrap'));
            wheel_move_iti_idx = find(ismember(wheel_starts, ...
                Timeline.rawDAQTimestamps(~stimOn_epochs)));

            % Get time from stim to rewarded movement onset and feedback
            stim_to_move = wheel_starts(wheel_move_stim_idx) - stimOn_times(1:n_trials);
            stim_to_feedback = signals_events.responseTimes(1:n_trials)' - stimOn_times(1:n_trials);

            % (choice and outcome)
            go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
            go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
            trial_choice = go_right(1:n_trials)' - go_left(1:n_trials)';
            trial_outcome = signals_events.hitValues(1:n_trials)' - signals_events.missValues(1:n_trials)';
            sides = [-1, 1];
            choices = [-1, 1];
            conditions = combvec(sides, choices)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.trialSideValues(1:n_trials)', ...
                trial_choice(1:n_trials)];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');

        case {'vanillaChoiceworld', 'vanillaChoiceworldBias', 'vanillaChoiceworldNoRepeats'}
            % Hit/miss recorded for last trial, circshift to align
            signals_events.hitValues = circshift(signals_events.hitValues, [0, -1]);
            signals_events.missValues = circshift(signals_events.missValues, [0, -1]);

            % Get number of completed trials (if uncompleted last trial)
            n_trials = length(signals_events.endTrialTimes);

            % Get stim on times by closest photodiode flip
            [~, closest_stimOn_photodiode] = ...
                arrayfun(@(x) min(abs(signals_events.stimOnTimes(x)- ...
                photodiode_flip_times)), ...
                1:n_trials);
            stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);

            % Check that the stim times aren't off by a certain threshold
            % (skip the first one - that's usually delayed a little)
            stim_time_offset_thresh = 0.1;
            if any(abs(stimOn_times(2:end)-signals_events.stimOnTimes(2:n_trials)') >= ...
                    stim_time_offset_thresh)
                figure;
                plot(stimOn_times-signals_events.stimOnTimes(1:n_trials)', '.k')
                line(xlim, repmat(stim_time_offset_thresh, 2, 1), 'color', 'r');
                line(xlim, repmat(-stim_time_offset_thresh, 2, 1), 'color', 'r');
                warning('Stim signals/photodiode offset over threshold');
                xlabel('Stim number');
                ylabel('Photodiode - signals stim time');
                title([animal, ' ', day, ' ', num2str(experiment)]);
            end

            % Get first movement time after stim onset
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);

            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)
            stim_to_move = padarray(wheel_move_time-stimOn_times, [n_trials - length(stimOn_times), 0], NaN, 'post');
            stim_to_feedback = signals_events.responseTimes(1:n_trials)' - stimOn_times(1:n_trials);

            % (early vs late move)
            trial_timing = 1 + (stim_to_move > 0.5);

            % (choice and outcome)
            go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
            go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
            trial_choice = go_right(1:n_trials)' - go_left(1:n_trials)';
            trial_outcome = signals_events.hitValues(1:n_trials)' - signals_events.missValues(1:n_trials)';

            % (trial conditions: [contrast,side,choice,timing])
            contrasts = [0, 0.06, 0.125, 0.25, 0.5, 1];
            sides = [-1, 1];
            choices = [-1, 1];
            timings = [1, 2];

            conditions = combvec(contrasts, sides, choices, timings)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.trialContrastValues(1:n_trials)', signals_events.trialSideValues(1:n_trials)', ...
                trial_choice(1:n_trials), trial_timing(1:n_trials)];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');
        case {'choiworldNoGoParameterHack_noWhiteNoise', 'noGo_stage4', 'noGo_stage4_q', 'noGo_stage5', 'noGo_stage6'} % stimType,
            % Hit/miss recorded for last trial, circshift to align
            response_trials = 1:length(block.events.endTrialValues);
            block.events.trialSideValues(response_trials) = 1;


            % Get number of completed trials (if uncompleted last trial)
            %keep pones with logged stimN (= not first and repeat on
            %incorrect)

            n_trials = [length(signals_events.stimulusOnTimes) - ...
                length(find(signals_events.stimulusOnTimes > signals_events.stimulusTypeTimes(1))), length(signals_events.endTrialTimes)];
            if n_trials(1) == 0
                n_trials = 1:n_trials(end);
            end
            % Get stim on times by closest photodiode flip
            [~, closest_stimOn_photodiode] = ...
                arrayfun(@(x) min(abs(signals_events.stimulusOnTimes(x)- ...
                photodiode_flip_times)), ...
                n_trials(1):n_trials(end));
            stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);


            % Check that the stim times aren't off by a certain threshold
            % (skip the first one - that's usually delayed a little)
            stim_time_offset_thresh = 0.1;
            if any(abs(stimOn_times(n_trials(2):n_trials(end))-signals_events.stimulusOnTimes(n_trials(2):n_trials(end))') >= ...
                    stim_time_offset_thresh)
                figure;
                plot(stimOn_times(n_trials)-signals_events.stimulusOnTimes(n_trials)', '.k')
                line(xlim, repmat(stim_time_offset_thresh, 2, 1), 'color', 'r');
                line(xlim, repmat(-stim_time_offset_thresh, 2, 1), 'color', 'r');
                warning('Stim signals/photodiode offset over threshold');
                xlabel('Stim number');
                ylabel('Photodiode - signals stim time');
                title([animal, ' ', day, ' ', num2str(experiment)]);
            end

            % Get first movement time after stim onset
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);

            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.2; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)
            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            stim_to_feedback = signals_events.responseTimes(n_trials(1):n_trials(end))' - stimOn_times;

            % (early vs late move)
            trial_timing = 1 + (stim_to_move > 0.5);

            % (choice and outcome)
            trial_choice = signals_events.responseValues(n_trials(1):n_trials(end));
            correctResp = nan(size(signals_events.stimulusTypeValues(n_trials(1):n_trials(end)), 2), 1);
            correctResp(signals_events.stimulusTypeValues(n_trials(1):n_trials(end)) == 1) = -1;
            correctResp(signals_events.stimulusTypeValues(n_trials(1):n_trials(end)) == 2) = -1;
            correctResp(signals_events.stimulusTypeValues(n_trials(1):n_trials(end)) == 3) = 0;
            trial_outcome =  block.events.responseValues(n_trials(1):n_trials(end))'== correctResp;

            imageN = unique(signals_events.stimulusTypeValues);
            choices = [-1, 0, 1];
            outcomes = [0, 1];
            conditions = combvec(imageN, choices, outcomes)';
            
            trial_conditions = ...
                [signals_events.stimulusTypeValues(n_trials(1):n_trials(end))', ...
                trial_choice', trial_outcome];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');
            stimIDs = signals_events.stimulusTypeValues;


        case {'vanillaChoiceworldImgs'}
            % Hit/miss recorded for last trial, circshift to align
            signals_events.hitValues = circshift(signals_events.hitValues, [0, -1]);
            signals_events.missValues = circshift(signals_events.missValues, [0, -1]);

            % Get number of completed trials (if uncompleted last trial)
            %keep pones with logged stimN (= not first and repeat on
            %incorrect)

            n_trials = [length(signals_events.stimOnTimes) - ...
                length(find(signals_events.stimOnTimes > signals_events.stimNTimes(1))),...
                length(signals_events.endTrialTimes)];

            % Get stim on times by closest photodiode flip
            [~, closest_stimOn_photodiode] = ...
                arrayfun(@(x) min(abs(signals_events.stimOnTimes(x)- ...
                photodiode_flip_times)), ...
                n_trials(1):n_trials(end));
            stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
            %stimOn
            % Check that the stim times aren't off by a certain threshold
            % (skip the first one - that's usually delayed a little)
            stim_time_offset_thresh = 0.1;
            if any(abs(stimOn_times(2:end)-signals_events.stimOnTimes(n_trials(1)+1:n_trials(2))') >= ...
                    stim_time_offset_thresh)
                figure;
                plot(stimOn_times(2:end)-signals_events.stimOnTimes(n_trials(1)+1:n_trials(2))', '.k')
                line(xlim, repmat(stim_time_offset_thresh, 2, 1), 'color', 'r');
                line(xlim, repmat(-stim_time_offset_thresh, 2, 1), 'color', 'r');
                warning('Stim signals/photodiode offset over threshold');
                xlabel('Stim number');
                ylabel('Photodiode - signals stim time');
                title([animal, ' ', day, ' ', num2str(experiment)]);
            end

            % Get first movement time after stim onset
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);

            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)
            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            stim_to_feedback = signals_events.responseTimes(n_trials(1):n_trials(end))' - stimOn_times;

            % (early vs late move)
            trial_timing = 1 + (stim_to_move > 0.5);

            % (choice and outcome)
            go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
            go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
            trial_choice = go_right(n_trials(1):n_trials(end))' - go_left(n_trials(1):n_trials(end))';
            trial_outcome = signals_events.hitValues(n_trials(1):n_trials(end))' - signals_events.missValues(n_trials(1):n_trials(end))';

            % (trial conditions: [contrast,side,choice,timing])
            %theseTrialsAnalyze = signals_events.stimNValues(1:n_trials)- stimOn_times; %for some reason, some of first aren't logged - drop them
            imageN = unique(signals_events.stimNValues);
            sides = [-1, 1];
            choices = [-1, 1];
            timings = [1, 2];
            %%hacky-need to chagnge in future and check - stimNMvalues probably end ones need to not be used =- because some first stimN values not logged
            %stimOn_times = stimOn_times(n_trials);
            %wheel_move_time = wheel_move_time(n_trials);
            % signals_events.responseTimes = signals_events.responseTimes(n_trials(1):n_trials(end));
            % trial_outcome = trial_outcome(n_trials(1):n_trials(end));

            conditions = combvec(imageN, sides, choices, timings)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.stimNValues(n_trials(1):n_trials(end))', signals_events.trialSideValues(n_trials(1):n_trials(end))', ...
                trial_choice, trial_timing];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');

        case {'locationWorld', 'noGo_stage1', 'noGo_stage2', 'noGo_stage3'}
            % Hit/miss recorded for last trial, circshift to align
            signals_events.hitValues = circshift(signals_events.hitValues, [0, -1]);
            signals_events.missValues = circshift(signals_events.missValues, [0, -1]);

            % Get number of completed trials (if uncompleted last trial)
            %keep pones with logged stimN (= not first and repeat on
            %incorrect)

            n_trials = [length(signals_events.stimOnTimes) - ...
                length(find(signals_events.stimOnTimes > signals_events.trialNumTimes(1))), length(signals_events.responseTimes)];
            n_trials(1) = n_trials(1) + 1;
            % Get stim on times by closest photodiode flip
            [~, closest_stimOn_photodiode] = ...
                arrayfun(@(x) min(abs(signals_events.stimOnTimes(x)- ...
                photodiode_flip_times)), ...
                n_trials(1):n_trials(end));
            stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);


            % Check that the stim times aren't off by a certain threshold
            % (skip the first one - that's usually delayed a little)
            stim_time_offset_thresh = 0.1;
            if any(abs(stimOn_times(2:end)-signals_events.stimOnTimes(n_trials(1)+1:n_trials(2))') >= ...
                    stim_time_offset_thresh)
                figure;
                plot(stimOn_times(2:end)-signals_events.stimOnTimes(n_trials(1)+1:n_trials(2))', '.k')
                line(xlim, repmat(stim_time_offset_thresh, 2, 1), 'color', 'r');
                line(xlim, repmat(-stim_time_offset_thresh, 2, 1), 'color', 'r');
                warning('Stim signals/photodiode offset over threshold');
                xlabel('Stim number');
                ylabel('Photodiode - signals stim time');
                title([animal, ' ', day, ' ', num2str(experiment)]);
            end

            % Get first movement time after stim onset
            surround_time = [-0.5, 5];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);

            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)
            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            stim_to_feedback = signals_events.responseTimes(n_trials(1):n_trials(end))' - stimOn_times;

            % (early vs late move)
            trial_timing = 1 + (stim_to_move > 0.5);

            % (choice and outcome)
            go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
            go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
                (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
            trial_choice = go_right(n_trials(1):n_trials(end))' - go_left(n_trials(1):n_trials(end))';
            trial_outcome = signals_events.hitValues(n_trials(1):n_trials(end))' - signals_events.missValues(n_trials(1):n_trials(end))';

            % (trial conditions: [contrast,side,choice,timing])
            %theseTrialsAnalyze = signals_events.stimNValues(1:n_trials)- stimOn_times; %for some reason, some of first aren't logged - drop them
            sides = [-1, 1];
            choices = [-1, 1];
            timings = [1, 2];
            %%hacky-need to chagnge in future and check - stimNMvalues probably end ones need to not be used =- because some first stimN values not logged
            %stimOn_times = stimOn_times(n_trials);
            %wheel_move_time = wheel_move_time(n_trials);
            % signals_events.responseTimes = signals_events.responseTimes(n_trials(1):n_trials(end));
            % trial_outcome = trial_outcome(n_trials(1):n_trials(end));

            [an, bn, cn] = ndgrid(sides, choices, timings); 
            conditions = [an(:), bn(:), cn(:)];

            %conditions = combvec(sides, choices, timings)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.trialSideValues(n_trials(1):n_trials(end))', ...
                trial_choice, trial_outcome];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');


      case 'AP_choiceWorldStimPassive'
            % This is kind of a dumb hack to get the stimOn times, maybe not
            % permanent unless it works fine: get stim times by checking for
            % close to the median photodiode flip difference
            block_stim_iti = mean(diff(block.stimWindowUpdateTimes));

            photodiode_flip_diff = diff(stimScreen_on_t(photodiode_flip));
            median_photodiode_flip_diff = mode(round(photodiode_flip_diff*10)/10);

            stimOn_idx = find(abs(photodiode_flip_diff-median_photodiode_flip_diff) < 0.1);

            stimOn_times = stimScreen_on_t(photodiode_flip(stimOn_idx))';

            % Set stimID as the contrast*side
            % (use last n values - sometimes short buffer times means some
            % stimuli in the beginning could be missed)
            use_signals_stim = size(signals_events.visualParamsValues, 2) - length(stimOn_times) + 1: ...
                size(signals_events.visualParamsValues, 2);
            stimIDs = sign(signals_events.visualParamsValues(1, use_signals_stim))' .* ...
                signals_events.visualParamsValues(2, use_signals_stim)';



        case {'JF_GratingPassive', 'JF_GratingPassiveVarITI', 'JF_GratingPassiveVarITI_moreComb', 'JF_GratingPassiveVarITI_moreCombnew', ...
                'JF_GratingPassiveVarITI_moreCombnew_correct'}
            stimOn_times = photodiode_flip_times(2:2:end);
            JF_correct_passive_photodiode; 

            % Get stim ID and conditions
            azimuths = unique(signals_events.stimAzimuthValues);
            spatialFreq = unique(signals_events.stimSpatialFreqValues);
            orientations = unique(signals_events.stimOrientationValues);

            conditions = combvec(azimuths, spatialFreq, orientations)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.stimAzimuthValues', signals_events.stimSpatialFreqValues', ...
                signals_events.stimOrientationValues'];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');
            [~, stimIDs] = ismember(trial_conditions, conditions, 'rows');

        case {'JF_Locations', 'JF_LocationsFit', 'JF_LocationsVarITI', 'JF_locations', 'JF_locationsFit', 'JF_locationsVarITI', ...
                'JF_locationsFitVarITIGrating'}
            stimOn_times = photodiode_flip_times(2:2:end);
            JF_correct_passive_photodiode; 

            % Get stim ID and conditions
            conditions = unique(signals_events.stim_idValues)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.stim_idValues'];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');
            [~, stimIDs] = ismember(trial_conditions, conditions, 'rows');

        case {'JF_natural_images', 'JF_natural_imagesVarITI', 'JF_natural_images_VarITInew', 'JF_natural_images_VarITI', 'JF_natural_imagesFitVarITI'}
            stimOn_times = photodiode_flip_times(2:2:end);
            JF_correct_passive_photodiode; 
            

            % Get stim ID and conditions
            conditions = unique(signals_events.stim_idValues)';
            n_conditions = size(conditions, 1);

            trial_conditions = ...
                [signals_events.stim_idValues'];
            [~, trial_id] = ismember(trial_conditions, conditions, 'rows');
            [~, stimIDs] = ismember(trial_conditions, conditions, 'rows');
            wheel_time = block.inputs.wheelTimes;
            wheel = block.inputs.wheelValues;
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            if size(stimOn_times, 2) > 1
                stimOn_times = permute(stimOn_times, [2, 1]);
            end
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);
            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.05;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;
            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            

        case {'JF_choiceworldStimuli', 'JF_choiceworldStimuli_wheel', 'JF_choiceworldStimuli_wheel_left_center', ...
                'JF_choiceworldStimuli_wheel_left_centerplok', 'JF_choiceworldStimuli_wheel_left_centerpl'}

            % Get stim times (first flip is initializing gray to black)
            stimOn_times = photodiode_flip_times(2:2:end);
            JF_correct_passive_photodiode; 

            if isfield(block.events, 'stim_aziValues')
                if max(signals_events.stim_idValues) > 44
                    conditions = unique([ceil(signals_events.stim_idValues/3)', block.events.stim_aziValues'], 'rows')';
                    n_conditions = size(conditions, 1);

                    trial_conditions = ... ,
                        [ceil(signals_events.stim_idValues/3)', block.events.stim_aziValues'];
                    [~, trial_id] = ismember(trial_conditions(:, 2), conditions(:, 2), 'rows');
                    [~, stimIDs] = ismember(trial_conditions(:, 1), conditions(:, 1), 'rows');

                else
                    conditions = unique([ceil(signals_events.stim_idValues)', block.events.stim_aziValues'], 'rows')';
                    n_conditions = size(conditions, 1);

                    trial_conditions = ... ,
                        [ceil(signals_events.stim_idValues)', block.events.stim_aziValues'];
                    trial_id = trial_conditions(:, 2);
                    stimIDs = trial_conditions(:, 1);
                end


            else
                warning('Azimuth not saved in the choiceworld stim version')
                conditions = unique(ceil(signals_events.stim_idValues/3))';
                tentativeAzi = zeros(size(signals_events.stim_idValues, 2), 1);
                tentativeAzi(signals_events.stim_idValues <= 22) = -90;
                tentativeAzi(signals_events.stim_idValues > 44) = 90;
                n_conditions = size(conditions, 1);
                trial_conditions = ...
                    [ceil(signals_events.stim_idValues/3)', tentativeAzi];
                [~, trial_id] = ismember(trial_conditions(:, 1), conditions, 'rows');
                [~, stimIDs] = ismember(trial_conditions(:, 1), conditions, 'rows');

            end


            wheel_time = block.inputs.wheelTimes;
            wheel = block.inputs.wheelValues;
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            if size(stimOn_times, 2) > 1
                stimOn_times = permute(stimOn_times, [2, 1]);
            end
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);
            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)

            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            no_move_trials = isnan(stim_to_move) | stim_to_move < 0.2 | stim_to_move > 0.2;

        case 'JF_choiceworldStimuli_wheel_left_center_all'


            stimOn_times = photodiode_flip_times(2:2:end);
            JF_correct_passive_photodiode; 

            % sanity check: times between stim on times in signals
            signals_photodiode_iti_diff = diff(signals_events.stimOnTimes(2:end)) - diff(stimOn_times) - 0.5';
            if any(signals_photodiode_iti_diff > 0.1)
                warning('mismatching signals/photodiode stim ITIs')
            end


            conditions = unique([signals_events.stim_idValues', block.events.stim_aziValues'], 'rows')';
            n_conditions = size(conditions, 1);

            trial_conditions = ... ,
                [ceil(signals_events.stim_idValues)', block.events.stim_aziValues'];
            trial_id = trial_conditions(:, 2);
            trial_conditions(trial_conditions(:, 1) > 13, 1) = trial_conditions(trial_conditions(:, 1) > 13, 1) - 13;
            stimIDs = trial_conditions(:, 1);


            wheel_time = block.inputs.wheelTimes;
            wheel = block.inputs.wheelValues;
            surround_time = [-0.5, 2];
            surround_sample_rate = 1 / Timeline.hw.samplingInterval; % (match this to framerate)
            surround_time_points = surround_time(1):1 / surround_sample_rate:surround_time(2);
            if size(stimOn_times, 2) > 1
                stimOn_times = permute(stimOn_times, [2, 1]);
            end
            pull_times = bsxfun(@plus, stimOn_times, surround_time_points);

            stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity, pull_times);
            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;

            % Get conditions for all trials

            % (trial_timing)

            stim_to_move = padarray(wheel_move_time-stimOn_times, [length(stimOn_times) - length(stimOn_times), 0], NaN, 'post');
            no_move_trials = isnan(stim_to_move) | stim_to_move < 0.2 | stim_to_move > 0.2;


        otherwise
            warning(['Signals protocol with no analysis script:', expDef]);
    end
