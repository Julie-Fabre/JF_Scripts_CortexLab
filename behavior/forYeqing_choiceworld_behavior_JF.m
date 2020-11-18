
%% Get and plot single mouse behavior (vanillaChoiceworld and choiceWorldJF)
corona = 0;
animalsAll = {'AP085', 'AP086', 'AP087'}; %chnage to the names of the animals you want here
for iAnimal = 1:size(animalsAll, 2)

    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'choiceworld'; %protocol name contains this name
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, protocol); %find the experiments with the choiceworld protocol

    bhv = struct; %initialize structure
 for curr_day = 1:length(experiments)

        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;

        % If multiple experiments, only use the last one (usually multiple
        % happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)

            experiment = experiment_num(curr_experiment);

            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');

            % Load the block file
            load(block_filename)
            if length(block.events.hitTimes) <=1
                removeThis(curr_day)=1;
            else
                removeThis(curr_day)=0;
            end
        end
 end
 experiments(logical(removeThis))=[];
    for curr_day = 1:length(experiments)

        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;

        % If multiple experiments, only use the last one (usually multiple
        % happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)

            experiment = experiment_num(curr_experiment);

            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');

            % Load the block file
            load(block_filename)
            if length(block.events.hitTimes) <=1
                experiments(curr_day)=[];
                curr_day = curr_day - 1;
            else

                % Get protocol name
                [~, protocol] = fileparts(block.expDef);

                % Time of session (in minutes)
                session_duration = block.duration / 60;

                % Trial counts
                n_trials = length(block.paramsValues);
                total_water = sum(block.outputs.rewardValues);

                % Estimate reaction time
                % (evenly resample velocity - not even natively);
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

                response_trials = 1:length(block.events.responseValues);
                if isfield(block.events, 'stimOnTimes')
                    trial_wheel_starts = arrayfun(@(x) ...
                        wheel_starts(find(wheel_starts > block.events.stimOnTimes(x), 1)), ...
                        response_trials);
                    trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
                    %  removeThis = find(cellfun(@isempty,trial_wheel_starts));
                else
                    % block.events.stimulusOnTimes = block.events.stimulusOnTimes(1:length(response_trials));
                    trial_wheel_starts = arrayfun(@(x) ...
                        wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                        response_trials, 'UniformOutput', false);
                    removeThis = find(cellfun(@isempty, trial_wheel_starts));
                    trial_wheel_starts = trial_wheel_starts(find(~cellfun(@isempty, trial_wheel_starts))); %remove empty, and back to vector
                    trial_wheel_starts = cell2mat(trial_wheel_starts);
                    trial_move_t = trial_wheel_starts - block.events.stimulusOnTimes(response_trials(1:end-length(removeThis)));
                end
                % Wheel movements/biases
                left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
                right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
                wheel_bias = (nansum(right_wheel_velocity) - nansum(left_wheel_velocity)) / ...
                    (nansum(right_wheel_velocity) + nansum(left_wheel_velocity));

                % Get reaction times for each stimulus
                if isfield(block.events, 'trialContrastValues')
                    trial_stim = block.events.trialContrastValues(response_trials) .* block.events.trialSideValues(response_trials);
                    stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1; 1], [], 1));
                    [~, trial_stim_idx] = ismember(trial_stim, stim_list);
                    stim_rxn_time = accumarray(trial_stim_idx(response_trials(1:end))', trial_move_t', [11, 1], @nanmedian, nan);

                else

                    trial_stim = block.events.contrastLeftValues(response_trials) - block.events.contrastRightValues(response_trials);
                    stim_list = unique(reshape(unique([block.events.contrastSetValues, -block.events.contrastSetValues]).*[-1; 1], [], 1));
                    [~, trial_stim_idx] = ismember(trial_stim, stim_list);
                    stim_rxn_time = accumarray(trial_stim_idx(response_trials(1:end-length(removeThis)))', trial_move_t', [11, 1], @nanmedian, nan);

                end


                % Performance (note that this excludes repeat on incorrect trials)
                if isfield(block.events, 'sessionPerformanceValues')
                    performance = block.events.sessionPerformanceValues(:, end-10:end);

                    %save performance
                    bhv.conditions = performance(1, :);
                    bhv.n_trials_condition(curr_day, :) = performance(2, :);
                    bhv.go_left_trials(curr_day, :) = performance(end, :);

                    if ~isempty(block.events.totalWaterValues)
                        bhv.rewardValue(curr_day) = block.events.totalWaterValues(1);
                    else
                        bhv.rewardValue(curr_day) = NaN;
                    end
                    bhv.stim_rxn_time(curr_day, :) = stim_rxn_time;
                else
                    % performance = block.events.sessionPerformanceValues(:,end-10:end);

                    %save performance
                    bhv.conditions = unique(reshape(unique([block.events.contrastSetValues, -block.events.contrastSetValues]).*[-1; 1], [], 1));
                    if numel(block.events.repeatNumValues) == numel(block.events.feedbackValues)
                        %repeatOnMisses = zeros(size(block.events.feedbackValues,2),1);
                        %repeatOnMisses = [block.events.repeatNumValues(2:end) > 1 & block.events.feedbackValues(1:end-1)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                    else
                        %repeatOnMisses = zeros(size(block.events.feedbackValues));all
                        %repeatOnMisses = [block.events.repeatNumValues(2:end-1) > 1 & block.events.feedbackValues(1:end-1)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                    end

                    for iCond = 1:size(bhv.conditions, 1)
                        bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond)));
                        bhv.n_trials_condition_not_repeat(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond))); %&~repeatOnMisses
                        bhv.go_left_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.responseValues == -1)); %&~repeatOnMisses
                    end
                    performance = bhv.go_left_trials ./ bhv.n_trials_condition_not_repeat; %removing repeat on miss trials

                    if contains(block.expDef, 'rew') %version of the task where manual rewards can be inputed
                        manualRew = length(strfind(block.inputs.keyboardValues, 'r'));
                    else
                        manualRew = 0;
                    end
                    rSize = block.paramsValues.rewardSize;
                    bhv.rewardValue(curr_day) = sum(block.outputs.rewardValues) - (manualRew * rSize);

                    bhv.stim_rxn_time(curr_day, :) = stim_rxn_time;


                end

                % Get whether all contrasts were used
                use_all_contrasts = all(block.events.useContrastsValues(end-10:end));

                % Store in behavior structure
                bhv.protocol{curr_day} = protocol;
                bhv.session_duration(curr_day) = session_duration;
                bhv.n_trials(curr_day) = n_trials;
                bhv.total_water(curr_day) = total_water;
                bhv.wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
                bhv.wheel_bias(curr_day) = wheel_bias;
                bhv.use_all_contrasts(curr_day) = use_all_contrasts;


                AP_print_progress_fractionJF(curr_day, length(experiments));
            end
        end
    end

    % Plot summary
    day_num = cellfun(@(x) datenum(x), {experiments.day});
    day_labels = cellfun(@(day, protocol) [protocol(19:end), ' ', day(6:end)], ...
        {experiments.day}, bhv.protocol, 'uni', false);
    figure('Name', animal)

    % Trials and water
    combine_days = bhv.use_all_contrasts;
    %combine_days=0;


    subplot(321)


    yyaxis left

    plot(day_num, bhv.n_trials, 'linewidth', 2);
    ylabel('Trials');
    yyaxis right

    plot(day_num, bhv.total_water, 'linewidth', 2);
    ax = gca;
    hold on;

    ylabel('Total water (ul)');
    xlabel('Session');
    set(gca, 'XTick', day_num);
    set(gca, 'XTickLabel', day_labels);
    set(gca, 'XTickLabelRotation', 90);

    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i), 1, 2), ylim, 'color', 'k');
    end

    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i), 1, 2), ylim, 'color', 'r', 'linestyle', '--');
    end

    % Wheel movement and bias

    subplot(322)


    yyaxis left

    plot(day_num, bhv.wheel_velocity, 'linewidth', 2);
    ylabel('Wheel movement');
    yyaxis right
    plot(day_num, bhv.wheel_bias, 'linewidth', 2);
    ylim([-1, 1]);
    line(xlim, [0, 0]);
    ylabel('Wheel bias');
    xlabel('Session');
    set(gca, 'XTick', day_num);
    set(gca, 'XTickLabel', day_labels);
    set(gca, 'XTickLabelRotation', 90);

    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i), 1, 2), ylim, 'color', 'k');
    end

    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i), 1, 2), ylim, 'color', 'r', 'linestyle', '--');
    end

    % Psychometric over all days

    subplot(323)


    im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials), bhv.go_left_trials./bhv.n_trials_condition, [0, 1]);
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    caxis = [0, 1];
    ylabel(c, 'Go left (frac)');
    xlabel('Condition');
    ylabel('Session');
    set(gca, 'YTick', 1:length(experiments));
    set(gca, 'YTickLabel', day_labels);
    axis square;
    hold on;
    if any([experiments.imaging])
        plot(0, find([experiments.imaging]), '.k');
    end
    if any([experiments.ephys])
        plot(0, find([experiments.ephys]), 'ok');
    end

    % Reaction time over days

    subplot(324)


    im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials), bhv.stim_rxn_time);
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    % caxis([0.2,0.8])
    c = colorbar;
    ylabel(c, 'Reaction time');
    xlabel('Condition');
    ylabel('Session');
    set(gca, 'YTick', 1:length(experiments));
    set(gca, 'YTickLabel', day_labels);
    axis square;
    hold on;
    if any([experiments.imaging])
        plot(0, find([experiments.imaging]), '.k');
    end
    if any([experiments.ephys])
        plot(0, find([experiments.ephys]), 'ok');
    end
    subplot(3, 2, 5);
    hold on;

    last_days_performance = bhv.go_left_trials(end, :) ./ bhv.n_trials_condition(end, :);
    errorbar(bhv.conditions, nanmean(last_days_performance, 1), ...
        nanstd(last_days_performance, [], 1)./sqrt(sum(~isnan(last_days_performance))), 'k', 'linewidth', 2);
    xlim([-1, 1]);
    ylim([0, 1]);
    line(xlim, [0.5, 0.5], 'color', 'k', 'linestyle', '--');
    line([0, 0], ylim, 'color', 'k', 'linestyle', '--');
    axis square;
    xlabel('Condition');
    ylabel('Fraction go left');

    % Reaction time of combined days that use all contrasts
    subplot(3, 2, 6);
    hold on;
    % combine_days = bhv.use_all_contrasts;
    last_days_rxn_time = bhv.stim_rxn_time(end, :);
    errorbar(bhv.conditions, nanmean(last_days_rxn_time, 1), ...
        nanstd(last_days_rxn_time, [], 1)./sqrt(sum(~isnan(last_days_rxn_time))), 'k', 'linewidth', 2);
    xlim([-1, 1]);
    line(xlim, [0.5, 0.5], 'color', 'k', 'linestyle', '--');
    axis square;
    xlabel('Condition');
    ylabel('Reaction time');
    % Psychometric of combined days that use all contrasts
    if sum(combine_days) ~= 0
        subplot(3, 2, 5);
        hold on;

        combine_days_performance = bhv.go_left_trials(combine_days, :) ./ bhv.n_trials_condition(combine_days, :);
        errorbar(bhv.conditions, nanmean(combine_days_performance, 1), ...
            nanstd(combine_days_performance, [], 1)./sqrt(sum(~isnan(combine_days_performance))), 'k', 'linewidth', 2);
        xlim([-1, 1]);
        ylim([0, 1]);
        line(xlim, [0.5, 0.5], 'color', 'k', 'linestyle', '--');
        line([0, 0], ylim, 'color', 'k', 'linestyle', '--');
        axis square;
        xlabel('Condition');
        ylabel('Fraction go left');

        % Reaction time of combined days that use all contrasts
        subplot(3, 2, 6);
        hold on;
        combine_days = bhv.use_all_contrasts;
        combine_days_rxn_time = bhv.stim_rxn_time(combine_days, :);
        errorbar(bhv.conditions, nanmean(combine_days_rxn_time, 1), ...
            nanstd(combine_days_rxn_time, [], 1)./sqrt(sum(~isnan(combine_days_rxn_time))), 'k', 'linewidth', 2);
        xlim([-1, 1]);
        line(xlim, [0.5, 0.5], 'color', 'k', 'linestyle', '--');
        axis square;
        xlabel('Condition');
        ylabel('Reaction time');
    end

end