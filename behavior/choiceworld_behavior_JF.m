
%% Get and plot single mouse behavior (vanillaChoiceworld and choiceWorldJF)

animalsAll = {'AL035'}; %JF029','JF032', 'JF033', 'JF034', 'JF035' };%{'JF032'};%{'JF025', 'JF026', 'JF028', 'JF029','JF032', 'JF033', 'JF034', 'JF035' };%chnage to the names of the animals you want here
for iAnimal = 1:size(animalsAll, 2)

    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'choiceworld'; %protocol name contains this name
    protocol2 = 'JF';
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, protocol2); %find the experiments with the choiceworld protocol
    %experiments(end)=[];
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
            %             figure();
            %             plot(wheel_values_resample(1:20000))
            %             figure();
            %             plot(wheel_velocity(1:20000))
            %            % wheel_t_resample %during all stims, average wheel movement
            %            if isfield(block.events, 'stimOnTimes')
            %                trial_wheel_movs_idx = arrayfun(@(x) ...
            %                     find(wheel_t_resample > block.events.stimOnTimes(x), 1), ...
            %                     response_trials);
            %            else
            %             trial_wheel_movs_idx = arrayfun(@(x) ...
            %                     find(wheel_t_resample > block.events.stimulusOnTimes(x), 1), ...
            %                     response_trials);
            %            end
            %            tt = nan(max(diff(trial_wheel_movs_idx)),size(trial_wheel_movs_idx,2)-1);
            %            trial_wheel_movs = nan(max(diff(trial_wheel_movs_idx)),size(trial_wheel_movs_idx,2)-1);
            %            for iW = 1:size(trial_wheel_movs_idx,2)-1
            %                tt(1:size(trial_wheel_movs_idx(iW):trial_wheel_movs_idx(iW+1),2),iW) = trial_wheel_movs_idx(iW):trial_wheel_movs_idx(iW+1);
            %                trial_wheel_movs(1:size(trial_wheel_movs_idx(iW):trial_wheel_movs_idx(iW+1),2),iW) = wheel_velocity(tt(1:size(trial_wheel_movs_idx(iW):trial_wheel_movs_idx(iW+1),2),iW));
            %            end
            %            figure();
            %            plot(wheel_t_resample(1:max(diff(trial_wheel_movs_idx)+1))-wheel_t_resample(1),nanmean(trial_wheel_movs'))
            %

            % padded array
            if isfield(block.paramsValues, 'noGoWheelGain')
                if ~isempty(find(block.events.orientationValues))
                    %                 nGwg = block.paramsValues.noGoWheelGain;
                    %                 gwg = block.paramsValues.wheelGain;
                    %                 %cell2mat( arrayfun(@(c) block.paramsValues.noGoWheelGain(2,:), (1:length(block.paramsValues)).', 'Uniform', 0) )
                    %                 go = block.events.orientationValues == -45 | block.events.orientationValues == 0;
                    %                 nogo =  block.events.orientationValues == 45 ;
                    %                 wheel_thresh = 0.005 .*(nogo.* nGwg) +...
                    %                     0.005 .* (go.* gwg) ;%QQ make dep on wheel gain
                    wheel_thresh = 0.025;


                end
            else
                if isfield(block.events, 'wheelGainValues')
                    wheel_thresh = 0.005 .* block.events.wheelGainValues(1);
                else
                    wheel_thresh = 0.025;
                end
            end
            %             figure();
            %             plot(wheel_velocity)
            %             hold on;
            %             scatter(trial_wheel_starts, ones(size(trial_wheel_starts,2),1))
            wheel_thresh = 0.025;
            wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
                abs(wheel_velocity(2:end)) > wheel_thresh);

            response_trials = 1:length(block.events.responseValues);
            if isfield(block.events, 'stimOnTimes')
                trial_wheel_starts = arrayfun(@(x) ...
                    wheel_starts(find(wheel_starts > block.events.stimOnTimes(x), 1)), ...
                    response_trials);
                trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
            elseif isfield(block.paramsValues, 'noGoWheelGain')
                try
                    trial_wheel_starts = arrayfun(@(x) ...
                        wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                        response_trials);
                catch
                    try
                        response_trials = response_trials(1:end-1);
                        trial_wheel_starts = arrayfun(@(x) ...
                            wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                            response_trials);
                    catch
                        try
                            response_trials = response_trials(1:end-1);
                            trial_wheel_starts = arrayfun(@(x) ...
                                wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                                response_trials);
                        catch
                            try
                                response_trials = response_trials(1:end-1);
                                trial_wheel_starts = arrayfun(@(x) ...
                                    wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                                    response_trials);
                            catch
                                try
                                    response_trials = response_trials(1:end-1);
                                    trial_wheel_starts = arrayfun(@(x) ...
                                        wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                                        response_trials);
                                catch
                                    continue
                                end
                            end
                        end
                    end
                end
                trial_move_t = trial_wheel_starts - block.events.stimulusOnTimes(response_trials);
            else

                trial_wheel_starts = arrayfun(@(x) ...
                    wheel_starts(find(wheel_starts > block.events.stimulusOnTimes(x), 1)), ...
                    response_trials);
                trial_move_t = trial_wheel_starts - block.events.stimulusOnTimes(response_trials);
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
            elseif isfield(block.events, 'contrastSetValues')
                trial_stim = block.events.contrastLeftValues(response_trials) - block.events.contrastRightValues(response_trials);
                stim_list = unique(reshape(unique([block.events.contrastSetValues, -block.events.contrastSetValues]).*[-1; 1], [], 1));


            elseif isfield(block.events, 'contrastRightValues') && (contains(block.expDef, 'one') || sum(block.events.contrastRightValues) == 0 || sum(block.events.contrastLeftValues) == 0)
                trial_stim = block.events.contrastLeftValues(response_trials) - block.events.contrastRightValues(response_trials);
                stim_list = nan(11, 1);
                if sum(block.events.contrastRightValues) == 0
                    s1 = unique(reshape(unique([block.events.contrastLeftValues]), [], 1));

                else
                    s1 = -unique(reshape(unique([block.events.contrastRightValues]), [], 1));

                end
                stim_list(1:size(s1, 1)) = s1(1:end);

            else
                trial_stim = block.events.contrastLeftValues(response_trials) - block.events.contrastRightValues(response_trials);
                stim_list = nan(11, 1);
                s1 = unique(reshape(unique([block.events.contrastLeftValues, -block.events.contrastLeftValues]).*[-1; 1], [], 1));
                stim_list(1:ceil((length(s1) - 1)/2)) = s1(1:ceil((length(s1) - 1)/2));
                stim_list(end-ceil((length(s1) - 1)/2)+1:end) = s1(end-ceil((length(s1) - 1)/2)+1:end);
                %stim_list((length(stim_list) - 1)/2+1) = s1(round((length(s1) - 1)/2+1));

            end


            [~, trial_stim_idx] = ismember(trial_stim, stim_list);
            stim_rxn_time = accumarray(trial_stim_idx(response_trials)', trial_move_t', [11, 1], @nanmedian, nan);


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
            elseif isfield(block.events, 'orientationValues') %different orientations
                orientations = unique(block.events.orientationValues);
                %save performance
                if isfield(block.events, 'contrastSetValues')
                    bhv.conditions = unique(reshape(unique([block.events.contrastSetValues, -block.events.contrastSetValues]).*[-1; 1], [], 1));
                else
                    bhv.conditions = stim_list;
                end
                %if numel(block.events.repeatNumValues) == numel(block.events.feedbackValues)
                %repeatOnMisses = zeros(size(block.events.feedbackValues,2),1);
                %repeatOnMisses = [block.events.repeatNumValues(2:end) > 1 & block.events.feedbackValues(1:end-1)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                %repeatOnMisses = zeros(size(block.events.feedbackValues));
                %repeatOnMisses = block.events.repeatNumValues((1:size(response_trials,2))) > 1; %& block.events.feedbackValues(1:end)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%

                %else
                %repeatOnMisses = zeros(size(block.events.feedbackValues));
                repeatOnMisses = block.events.repeatNumValues(response_trials) > 1; %& block.events.feedbackValues(1:end)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                %end

                for iCond = 1:size(bhv.conditions, 1)
                    bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond)));
                    %                     if numel(block.events.orientationValues) == numel(trial_stim)
                    %                         bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues == 0 & ~repeatOnMisses));
                    %                         bhv.go_left_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == 0 & block.events.responseValues == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.n_trials_condition45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues == -45 & ~repeatOnMisses));
                    %                         bhv.go_left_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == -45 & block.events.responseValues == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.n_trials_conditionNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues == 45 & ~repeatOnMisses));
                    %                         bhv.go_left_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == 45 & block.events.responseValues == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == 0 & block.events.responseValues == 0)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == -45 & block.events.responseValues == 0)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues == 45 & block.events.responseValues == 0)); %&~repeatOnMisses
                    %
                    %                     elseif numel(trial_stim)+1 == numel(repeatOnMisses)
                    bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(response_trials) == 0 & ~repeatOnMisses));
                    bhv.go_left_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & ...
                        block.events.orientationValues(response_trials) == 0 & block.events.responseValues(response_trials) == -1)); %&~repeatOnMisses

                    bhv.n_trials_condition45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(response_trials) == -45 & ~repeatOnMisses));
                    bhv.go_left_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(response_trials) == -45 & block.events.responseValues(response_trials) == -1)); %&~repeatOnMisses

                    bhv.n_trials_conditionNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(response_trials) == 45 & ~repeatOnMisses));
                    bhv.go_left_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(response_trials) == 45 & block.events.responseValues(response_trials) == -1)); %&~repeatOnMisses
                    %rt = trial_move_t(response_trials);
                    %min(block.events.responseTimes-block.events.contrastLeftTimes(1:end-1))
                    bb = block.paramsValues.noGoResponseWindow;
                    bhv.no_go_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(response_trials) == 0 ...
                        & (block.events.responseValues(response_trials) == 0 | block.events.responseTimes(response_trials) - block.events.contrastLeftTimes(response_trials) - 0.4 > bb))); %&~repeatOnMisses

                    bhv.no_go_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(response_trials) == -45 ...
                        & (block.events.responseValues(response_trials) == 0 | block.events.responseTimes(response_trials) - block.events.contrastLeftTimes(response_trials) - 0.4 > bb))); %&~repeatOnMisses

                    bhv.no_go_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(response_trials) == 45 ...
                        & (block.events.responseValues(response_trials) == 0 | block.events.responseTimes(response_trials) - block.events.contrastLeftTimes(response_trials) - 0.4 > bb))); %&~repeatOnMisses


                    %                     else
                    %                         bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(1:end-1) == 0 & ~repeatOnMisses));
                    %                         bhv.go_left_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == 0 & block.events.responseValues(1:size(response_trials,2)) == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.n_trials_condition45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(1:end-1) == -45 & ~repeatOnMisses));
                    %                         bhv.go_left_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == -45 & block.events.responseValues(1:size(response_trials,2)) == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.n_trials_conditionNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & block.events.orientationValues(1:end-1) == 45 & ~repeatOnMisses));
                    %                         bhv.go_left_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == 45 & block.events.responseValues(1:size(response_trials,2)) == -1)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == 0 & block.events.responseValues(1:size(response_trials,2)) == 0)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trials45(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == -45 & block.events.responseValues(1:size(response_trials,2)) == 0)); %&~repeatOnMisses
                    %
                    %                         bhv.no_go_trialsNoGo(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.orientationValues(1:end-1) == 45 & block.events.responseValues(1:size(response_trials,2)) == 0)); %&~repeatOnMisses
                    %
                    %end
                end
                performance = bhv.go_left_trials ./ bhv.n_trials_condition; %removing repeat on miss trials
                performance45 = bhv.go_left_trials45 ./ bhv.n_trials_condition45;
                performanceNoGo = bhv.go_left_trialsNoGo ./ bhv.n_trials_conditionNoGo;
                if contains(block.expDef, 'rew') %version of the task where manual rewards can be inputed
                    manualRew = length(strfind(block.inputs.keyboardValues, 'r'));
                else
                    manualRew = 0;
                end
                rSize = block.paramsValues.rewardSize;
                bhv.rewardValue(curr_day) = sum(block.outputs.rewardValues) - (manualRew * rSize);

                %                 if numel(block.events.orientationValues) == numel(trial_stim)
                %                     stim_rxn_time = accumarray(trial_stim_idx(block.events.orientationValues(1:end) == 0 & ~repeatOnMisses)',...
                %                         trial_move_t(block.events.orientationValues(1:end) == 0 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);
                %                     stim_rxn_time45 = accumarray(trial_stim_idx(block.events.orientationValues(1:end) == -45 & ~repeatOnMisses)',...
                %                         trial_move_t(block.events.orientationValues(1:end) == -45 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);
                %                     stim_rxn_timeNoGo = accumarray(trial_stim_idx(block.events.orientationValues(1:end) == 45 & ~repeatOnMisses)',...
                %                         trial_move_t(block.events.orientationValues(1:end) == 45 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);
                %
                %                 else
                stim_rxn_time = accumarray(trial_stim_idx(block.events.orientationValues(response_trials) == 0 & ~repeatOnMisses)', ...
                    trial_move_t(block.events.orientationValues(response_trials) == 0 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);
                stim_rxn_time45 = accumarray(trial_stim_idx(block.events.orientationValues(response_trials) == -45 & ~repeatOnMisses)', ...
                    trial_move_t(block.events.orientationValues(response_trials) == -45 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);
                stim_rxn_timeNoGo = accumarray(trial_stim_idx(block.events.orientationValues(response_trials) == 45 & ~repeatOnMisses)', ...
                    trial_move_t(block.events.orientationValues(response_trials) == 45 & ~repeatOnMisses)', [11, 1], @nanmedian, nan);

                %trial_stim_idx(response_trials)
                %                 end
                bhv.stim_rxn_time(curr_day, :) = stim_rxn_time;
                bhv.stim_rxn_time45(curr_day, :) = stim_rxn_time45;
                bhv.stim_rxn_timeNoGo(curr_day, :) = stim_rxn_timeNoGo;


            else
                % performance = block.events.sessionPerformanceValues(:,end-10:end);

                %save performance
                if isfield(block.events, 'contrastSetValues')
                    bhv.conditions = unique(reshape(unique([block.events.contrastSetValues, -block.events.contrastSetValues]).*[-1; 1], [], 1));
                else
                    bhv.conditions = stim_list;
                end
                if numel(block.events.repeatNumValues) == numel(block.events.feedbackValues)
                    %repeatOnMisses = zeros(size(block.events.feedbackValues,2),1);
                    %repeatOnMisses = [block.events.repeatNumValues(2:end) > 1 & block.events.feedbackValues(1:end-1)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                    repeatOnMisses = zeros(size(block.events.feedbackValues));
                    repeatOnMisses = block.events.repeatNumValues(1:end) > 1; %& block.events.feedbackValues(1:end)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%

                else
                    repeatOnMisses = zeros(size(block.events.feedbackValues));
                    repeatOnMisses = block.events.repeatNumValues(1:end-1) > 1; %& block.events.feedbackValues(1:end)==-1, 0];%zeros(size(block.events.feedbackValues,1),1);%
                end

                for iCond = 1:size(bhv.conditions, 1)
                    %bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond)));
                    bhv.n_trials_condition(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses));
                    bhv.go_left_trials(curr_day, iCond) = numel(find(trial_stim == bhv.conditions(iCond) & ~repeatOnMisses & block.events.responseValues == -1)); %&~repeatOnMisses
                end
                performance = bhv.go_left_trials ./ bhv.n_trials_condition; %removing repeat on miss trials

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
            if isfield(block.events, 'useContrastValues')
                use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            else
                use_all_contrasts = sum(~isnan(stim_list)) == 11;
            end
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

    inver = strcmp(bhv.protocol, 'oneSideChoiceWorldJF');
    % Plot summary
    day_num = cellfun(@(x) datenum(x), {experiments.day});
    day_labels = cellfun(@(day) [' ', day(6:end)], ...
        {experiments.day}, 'uni', false);
    figure('Name', animal)

    % Trials and water
    combine_days = bhv.use_all_contrasts;
    %combine_days=0;

    if isfield(bhv, 'go_left_trials45')
        subplot(431)
    else
        subplot(321)
    end

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

    if isfield(bhv, 'go_left_trials45')
        subplot(433)
    else
        subplot(322)
    end


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

    if isfield(bhv, 'go_left_trials45')
        subplot(434)
    else
        subplot(323)
    end

    bhv.go_left_trials(:, 3:9) = 0;
    bhv.n_trials_condition(:, 3:9) = 0;
    bhv.conditions(isnan(bhv.conditions)) = 0;
    im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials), abs(repmat(inver', 1, 11)-bhv.go_left_trials./bhv.n_trials_condition));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'Go left (frac)');
    xlabel('Condition');
    ylabel('Session');
    set(gca, 'YTick', 1:length(experiments));
    set(gca, 'YTickLabel', day_labels);
    axis square;
    caxis([0, 1])
    hold on;
    if any([experiments.imaging])
        plot(0, find([experiments.imaging]), '.k');
    end
    if any([experiments.ephys])
        plot(0, find([experiments.ephys]), 'ok');
    end

    % Reaction time over days

    if isfield(bhv, 'go_left_trials45')
        subplot(436)
    else
        subplot(324)
    end

    bhv.stim_rxn_time(:, 3:9) = NaN;
    % bhv.n_trials_condition(:,4:9) = 0;
    im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials), bhv.stim_rxn_time);
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    %caxis([0, 1])
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

    bhv.no_go_trials(:, 3:9) = 0;
    bhv.n_trials_condition(:, 3:9) = 0;
    if isfield(bhv, 'go_left_trials45')
        subplot(435)
        bhv.conditions(isnan(bhv.conditions)) = 0;
        im = imagesc(bhv.conditions, 1:size(bhv.no_go_trials), bhv.no_go_trials./bhv.n_trials_condition);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        c = colorbar;
        ylabel(c, 'No go (frac)');
        xlabel('Condition');
        ylabel('Session');
        set(gca, 'YTick', 1:length(experiments));
        set(gca, 'YTickLabel', day_labels);
        axis square;
        caxis([0, 1])
        hold on;
        if any([experiments.imaging])
            plot(0, find([experiments.imaging]), '.k');
        end
        if any([experiments.ephys])
            plot(0, find([experiments.ephys]), 'ok');
        end


        subplot(434)
        title('Go1')

        subplot(437)
        bhv.conditions(isnan(bhv.conditions)) = 0;
        im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials45), bhv.go_left_trials45./bhv.n_trials_condition45);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        c = colorbar;
        ylabel(c, 'Go left (frac)');
        xlabel('Condition');
        ylabel('Session');
        set(gca, 'YTick', 1:length(experiments));
        set(gca, 'YTickLabel', day_labels);
        axis square;
        caxis([0, 1])
        hold on;
        if any([experiments.imaging])
            plot(0, find([experiments.imaging]), '.k');
        end
        if any([experiments.ephys])
            plot(0, find([experiments.ephys]), 'ok');
        end
        title('Go2')

        subplot(438)
        bhv.conditions(isnan(bhv.conditions)) = 0;
        im = imagesc(bhv.conditions, 1:size(bhv.no_go_trials45), bhv.no_go_trials45./bhv.n_trials_condition45);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        c = colorbar;
        ylabel(c, 'No go (frac)');
        xlabel('Condition');
        ylabel('Session');
        set(gca, 'YTick', 1:length(experiments));
        set(gca, 'YTickLabel', day_labels);
        axis square;
        caxis([0, 1])
        hold on;
        if any([experiments.imaging])
            plot(0, find([experiments.imaging]), '.k');
        end
        if any([experiments.ephys])
            plot(0, find([experiments.ephys]), 'ok');
        end


        subplot(439)
        bhv.stim_rxn_time45(~any(bhv.stim_rxn_time45, 2), :) = NaN;
        im = imagesc(bhv.conditions, 1:size(bhv.go_left_trials45), bhv.stim_rxn_time45);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        %caxis([0, 01])
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


        subplot(4, 3, 10)
        bhv.conditions(isnan(bhv.conditions)) = 0;
        im = imagesc(bhv.conditions, 1:size(bhv.go_left_trialsNoGo), bhv.go_left_trialsNoGo./bhv.n_trials_conditionNoGo);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        c = colorbar;
        ylabel(c, 'Go left (frac)');
        xlabel('Condition');
        ylabel('Session');
        set(gca, 'YTick', 1:length(experiments));
        set(gca, 'YTickLabel', day_labels);
        axis square;
        %caxis([0, 1])
        hold on;
        if any([experiments.imaging])
            plot(0, find([experiments.imaging]), '.k');
        end
        if any([experiments.ephys])
            plot(0, find([experiments.ephys]), 'ok');
        end
        title('NoGo')

        subplot(4, 3, 11)
        bhv.conditions(isnan(bhv.conditions)) = 0;
        im = imagesc(bhv.conditions, 1:size(bhv.no_go_trialsNoGo), bhv.no_go_trialsNoGo./bhv.n_trials_conditionNoGo);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        c = colorbar;
        ylabel(c, 'No go (frac)');
        xlabel('Condition');
        ylabel('Session');
        set(gca, 'YTick', 1:length(experiments));
        set(gca, 'YTickLabel', day_labels);
        axis square;
        caxis([0, 1])
        hold on;
        if any([experiments.imaging])
            plot(0, find([experiments.imaging]), '.k');
        end
        if any([experiments.ephys])
            plot(0, find([experiments.ephys]), 'ok');
        end

        subplot(4, 3, 12)
        bhv.stim_rxn_timeNoGo(~any(bhv.stim_rxn_timeNoGo, 2), :) = NaN;
        im = imagesc(bhv.conditions, 1:size(bhv.go_left_trialsNoGo), bhv.stim_rxn_timeNoGo);
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);
        colormap(brewermap([], '*RdBu'));
        caxis([0, 1])
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


    else

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
end