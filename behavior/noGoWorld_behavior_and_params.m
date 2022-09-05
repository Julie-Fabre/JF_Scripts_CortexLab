% get parameters and plot: gain, no go thresh, go thresh
% p.dispThresholdNoGo = 5;
% p.wheelGain = 4;
% p.noGoWheelGain = 4;
% diff go1/no go no go percentage
% p.interactiveDelay = 0;
% p.preStimulusDelay = [0.5, 0.5, 0.5]';


function bhvOut = noGoWorld_behavior_and_params(animalsAll)

bhvOut = struct;
for iAnimal = 1:size(animalsAll, 2)


    animal = animalsAll{1, iAnimal}; %this animal
    if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
        thisSide = -1;
    else
        thisSide = 1;
    end
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
    theseD = 1:length(experiments);

    for curr_day = theseD


        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        trNum = [];
        % If multiple experiments, only use the largest one (usually multiple
        % happens if mess ups first/last one is good)
        %if curr_day ~= length(experiments)
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
        end

        %
        for curr_experiment = find(trNum == max(trNum))
            try
                experiment = experiment_num(curr_experiment);

                [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
                load(block_filename)

                [param_filename, param_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'parameters');
                load(param_filename)
                if isfield(parameters, 'dispThresholdNoGo')
                    bhv.dispThresholdNoGo(curr_day) = parameters.dispThresholdNoGo;
                    bhv.interactiveDelay(curr_day) = parameters.interactiveDelay;
                    bhv.noGoWheelGain(curr_day) = parameters.noGoWheelGain;
                else
                    bhv.dispThresholdNoGo(curr_day) = NaN;
                    bhv.interactiveDelay(curr_day) = NaN;
                    bhv.noGoWheelGain(curr_day) = NaN;
                end
                if isfield(parameters, 'diffWheelMove')
                    bhv.diffWheelMove(curr_day) = parameters.diffWheelMove;
                else
                    bhv.diffWheelMove(curr_day) = NaN;
                end
                if isfield(parameters, 'preStimulusDelay')
                    bhv.preStimulusDelay(curr_day) = unique(paramaters.preStimulusDelay);
                else
                    bhv.preStimulusDelay(curr_day) = NaN;
                end
                
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


                bhv.nTrials(curr_day, :) = nTrials;
                bhv.n_trials(curr_day, :) = n_trials;
                bhv.goLeft(curr_day, :) = goLeft;
                bhv.noGo(curr_day, :) = noGo;

                keep_day = [keep_day, curr_day];
                expDefNameStart = strfind(block.expDef, 'stage');
                expDefName = block.expDef(expDefNameStart:expDefNameStart+5);
                bhv.expDefName{curr_day} = expDefName;

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

%% summary: % no go per stim per mouse on last day only
plot()

end