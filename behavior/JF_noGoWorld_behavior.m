%clear all;
%clc;
function bhvOut = JF_noGoWorld_behavior(animalsAll)
plotAll = false;
animalsPhase1 = {[1, 2], [1, 2], [1:4]};
animalsPhase2 = {[3, 4], [3:5], [5:9]};
animalsPhase3 = {[5], [6], [10]};
animalsPhase4 = {[6:12], [7:14], [11:13]};
animalsPhase5 = {[13, 14, 16], [16], []};
%animalsAll = {'JF067'};
bhvOut = struct;
for iAnimal = 1:size(animalsAll, 2)
%bhv
myPaths;
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
    if strcmp(animalsAll{iAnimal}, 'JF067')
        theseD(26) = [];
    end
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
                
                go1TrialsFirstHalf = go1Trials;
                go1TrialsFirstHalf(round(length(go1Trials)/2)+1:end) = 0; 
                go2TrialsFirstHalf = go2Trials;
                go2TrialsFirstHalf(round(length(go1Trials)/2)+1:end) = 0; 
                noGoTrialsFirstHalf = noGoTrials;
                noGoTrialsFirstHalf(round(length(go1Trials)/2)+1:end) = 0; 
                
                go1TrialsSecHalf = go1Trials;
                go1TrialsSecHalf(1:round(length(go1Trials)/2)) = 0; 
                go2TrialsSecHalf = go2Trials;
                go2TrialsSecHalf(1:round(length(go1Trials)/2)) = 0; 
                noGoTrialsSecHalf = noGoTrials;
                noGoTrialsSecHalf(1:round(length(go1Trials)/2)) = 0;
                
                if contains(block.expDef, 'Hack') || contains(block.expDef, 'noGo_stage4') || contains(block.expDef, 'noGo_stage5')
                    if strcmp(animal, 'JF042') || strcmp(animal, 'JF043') || strcmp(animal, 'JF044')
                        val = 1;
                    else
                    val = -1 ;
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
                    
                nTrials1((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                    go1TrialsFirstHalf ...
                    );
                nTrials1((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                    go2TrialsFirstHalf ...
                    );
                nTrials1((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                    noGoTrialsFirstHalf ...
                    );
                nTrials2((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                    go1TrialsSecHalf ...
                    );
                nTrials2((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                    go2TrialsSecHalf ...
                    );
                nTrials2((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                    noGoTrialsSecHalf ...
                    );
                
                    
                    goLeft1((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go1TrialsFirstHalf ...
                        );
                    goLeft1((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go2TrialsFirstHalf ...
                        );
                    goLeft1((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & noGoTrialsFirstHalf ...
                        );
                    goLeft2((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go1TrialsSecHalf ...
                        );
                    goLeft2((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & go2TrialsSecHalf ...
                        );
                    goLeft2((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                        ...
                        block.events.hitValues(response_trials) == val & noGoTrialsSecHalf ...
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
                
                
                
                noGo1((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & go1TrialsFirstHalf ...
                    );
                noGo1((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & go2TrialsFirstHalf ...
                    );
                noGo1((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & noGoTrialsFirstHalf ...
                    );
                
                noGo2((iImg - 1)*3+1) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & go1TrialsSecHalf ...
                    );
                noGo2((iImg - 1)*3+2) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & go2TrialsSecHalf ...
                    );
                noGo2((iImg - 1)*3+3) = sum(~repeatOnMisses(response_trials) & ...
                    ...
                    (block.events.hitValues(response_trials) == 0) & noGoTrialsSecHalf ...
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
            surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
            if isfield(block.events, 'stimOnTimes')
                pull_times = bsxfun(@plus, block.events.stimOnTimes', surround_time_points);
            else
                pull_times = bsxfun(@plus, block.events.stimulusOnTimes', surround_time_points);
            end

            stim_aligned_wheel = interp1(wheel_t_resample, ...
                wheel_velocity, pull_times);
            
           
%     wheel_values_resample(wheel_values_resample> 2^31) = wheel_values_resample(wheel_values_resample> 2^31) - 2^32;
%     
%     [wheel_velocity,wheel_move,wheel_velocity_split] = ...
%     AP_parse_wheel(wheel_values_resample,wheel_resample_rate);
% 
%             wheel_click2mm = 0.4869; % lilrig's encoder
%         wheel_mm2deg = 4;  % gain in this expDef
%         wheel_position_mm = wheel_values_resample*wheel_click2mm;
%         wheel_position_deg = wheel_position_mm*wheel_mm2deg;
%         
%             stim_aligned_wheel_deg = interp1(wheel_t_resample, ...
%                 wheel_position_deg, pull_times);

            % (set a threshold in speed and time for wheel movement)
            thresh_displacement = 0.025;%% QQ SHOULD BE 0.025
            time_over_thresh = 0.05; % ms over velocity threshold to count
            samples_over_thresh = time_over_thresh .* surround_sample_rate;
            wheel_over_thresh_fullconv = convn( ...
                abs(stim_aligned_wheel) > thresh_displacement, ...
                ones(1, samples_over_thresh)) >= samples_over_thresh;
            wheel_over_thresh = wheel_over_thresh_fullconv(:, end-size(stim_aligned_wheel, 2)+1:end);

            [move_trial, wheel_move_sample] = max(wheel_over_thresh, [], 2);
            wheel_move_time = arrayfun(@(x) pull_times(x, wheel_move_sample(x)), 1:size(pull_times, 1))';
            wheel_move_time(~move_trial) = NaN;
            
            
            

                

                %trial_stim = block.events.trialContrastValues(response_trials).*block.events.trialSideValues(response_trials);
                %stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1;1],[],1));
                %[~,trial_stim_idx] = ismember(trial_stim,stim_list);
                if isfield(block.events, 'stimOnTimes')
                    stim_to_move = padarray(wheel_move_time'-block.events.stimOnTimes, ...
                        [n_trials - length(block.events.stimOnTimes), 0], NaN, 'post');
            
                    trial_wheel_starts = arrayfun(@(x) ...
                        wheel_starts(find(wheel_starts > block.events.stimOnTimes(x), 1)), ...
                        response_trials);

                    trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
                    stim_rxn_time = nanmedian(trial_move_t);
                    stim_rxn_timeSEM = nanstd(trial_move_t) ./ sqrt(numel(trial_move_t));

%                     stim_leeway = 0.1;
%             wheel_move_alt_stim_idx = ...
%                 arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
%                 cell2mat(alt_stimOn_times));
%             
%             alt_stim_to_move = ...
%                 mat2cell(wheel_starts(wheel_move_alt_stim_idx) - cell2mat(alt_stimOn_times), ...
%                 cellfun(@length,alt_stimOn_times));

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
                    bhv.stim_aligned_wheelGo1{curr_day} = stim_aligned_wheel(go1Trials,:);
                    bhv.stim_aligned_wheelGo1_inco{curr_day} = stim_aligned_wheel(go1Trials & block.events.hitValues(response_trials) ~= val,:);
                    bhv.stim_aligned_wheelGo1_corr{curr_day} = stim_aligned_wheel(go1Trials & block.events.hitValues(response_trials) == val,:);
                else
                    bhv.stim_aligned_wheelGo1{curr_day} = stim_aligned_wheel;
                    bhv.stim_aligned_wheelGo1_corr{curr_day} = [];
                    bhv.stim_aligned_wheelGo1_inco{curr_day} = [];
                end
                if sum(go2Trials) > 0
                    bhv.stim_aligned_wheelGo2{curr_day} = stim_aligned_wheel(go2Trials,:);
                    bhv.stim_aligned_wheelGo2_inco{curr_day} = stim_aligned_wheel(go2Trials & block.events.hitValues(response_trials) ~= val,:);
                    bhv.stim_aligned_wheelGo2_corr{curr_day} = stim_aligned_wheel(go2Trials & block.events.hitValues(response_trials) == val,:);
                
                else
                    bhv.stim_aligned_wheelGo2{curr_day} = [];
                    bhv.stim_aligned_wheelGo2_corr{curr_day} = [];
                    bhv.stim_aligned_wheelGo2_inco{curr_day} = [];
                end
                if sum(noGoTrials) > 0 
                    bhv.stim_aligned_wheelNoGo{curr_day} = stim_aligned_wheel(noGoTrials,:);
                    bhv.stim_aligned_wheelNoGo_corr{curr_day} = stim_aligned_wheel(noGoTrials & block.events.hitValues(response_trials) == 0,:);
                    bhv.stim_aligned_wheelNoGo_inco{curr_day} = stim_aligned_wheel(noGoTrials & block.events.hitValues(response_trials) ~= 0,:);
              
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
                if size(stim_to_move,1)==2
                    stim_to_move = stim_to_move(1,:);
                end
                bhv.stim_to_move{curr_day,1} = stim_to_move; 
                bhv.stim_to_moveMean(curr_day,1) = nanmean(stim_to_move);
                bhv.stim_to_moveMean(curr_day,2) = NaN;
                bhv.stim_to_moveMean(curr_day,3) = NaN;
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
                    bhv.stim_to_moveMean(curr_day,1) = nanmean(stim_to_move(go1Trials(1:end))); 
                    if any(go2Trials)
                        bhv.stim_to_move{curr_day,2} = stim_to_move(go2Trials(1:end) );
                         bhv.stim_to_move{curr_day,3} = stim_to_move(noGoTrials(1:end) ); 
                    elseif any(noGoTrials)
                         bhv.stim_to_move{curr_day,2} =  nan(size( movingFrac,2),1)';
                         bhv.stim_to_move{curr_day,3} = stim_to_move(noGoTrials(1:end) ); 
                    else
                        bhv.stim_to_move{curr_day,2} =  nan(size( movingFrac,2),1)';
                        bhv.stim_to_move{curr_day,3} = nan(size( movingFrac,2),1)';
                    end
                    bhv.stim_to_moveMean(curr_day,3) = nanmean(stim_to_move(noGoTrials(1:end))); 
                    
                    bhv.stim_to_move{curr_day,1} = stim_to_move(go1Trials(1:end) ); 
                   % bhv.stim_to_move{curr_day,2} = stim_to_move(go2Trials(1:end) ); 
                   % bhv.stim_to_move{curr_day,3} = stim_to_move(noGoTrials(1:end) ); 
                    
                    
                    bhv.movingFracGo1(curr_day, :) = movingFracGo1;
                    bhv.movingFracNoGo(curr_day, :) = movingFracNoGo;
                    bhv.trialMoveNoGoCorrect{curr_day, :, :} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                    bhv.trialMoveNoGoIncorrect{curr_day, :, :} = trial_wheel_move(noGoTrials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                    bhv.trialMoveGo1Correct{curr_day, :, :} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                    bhv.trialMoveGo2Correct{curr_day, :, :} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 1), :);
                    bhv.trialMoveGo1Incorrect{curr_day, :, :} = trial_wheel_move(go1Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                    bhv.trialMoveGo2Incorrect{curr_day, :, :} = trial_wheel_move(go2Trials(1:end) & (block.events.hitValues(response_trials(1:end)) == 0), :);
                else
                    %bhv.movingFrac(curr_day, 1) = movingFrac;
                    bhv.stim_to_move{curr_day,1} = stim_to_move;
                    bhv.stim_to_move{curr_day,2} = nan(size( movingFrac,2),1)';
                    bhv.stim_to_move{curr_day,3} = nan(size( movingFrac,2),1)';


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
        {experiments(keep_day).day},bhv.expDefName(keep_day), 'uni', false);
    bhvOut(iAnimal). dates =  cellfun(@(day, protocol) [day], ...
        {experiments(keep_day).day},bhv.expDefName(keep_day), 'uni', false);
    %  day_labels = cellfun(@(day, protocol) [protocol, ' ', day(6:end)], ...
    %    {experiments(keep_day).day},bhv.expDefName(keep_day), 'uni', false);
    %
    tempLabel = "\\color[rgb]{%s}";
       
 

    [unique_protocols,~,protocol_idx] = unique(bhv.expDefName(keep_day));
    protocol_col = hsv(length(unique_protocols));
    for iDay = 1:length(protocol_idx)
        thisColor = sprintf(tempLabel,num2str(protocol_col(protocol_idx(iDay),:)));
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
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',45);
    %set(gca, 'Color', 
    makepretty;

% if ismember(iAnimal, 1:3)
%     bhv.goLeft(keep_day, 3) = bhv.goLeft(keep_day, 3) -0.2;
%     bhv.noGo(keep_day, 3) = bhv.noGo(keep_day, 3) +0.2;
% end
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
    for iDay = 1:size(bhv.stim_to_move,1)
    [h1(iDay,1,:), bins1(iDay,1,:) ]= hist(bhv.stim_to_move{iDay,1}, 0:0.05:1.8);
    %day_idx_go1=[day_idx_go1; ones(size(bhv.stim_to_move{iDay,1},2),1)*iDay];
    n_rxn_altsample = 1000;
    rxn_measured_med(iDay,1) = nanmedian(bhv.stim_to_move{iDay,1});
    rxn_alt_med(iDay,1) = nanmedian(datasample(bhv.stim_to_move{iDay,1},n_rxn_altsample));
    [h1(iDay,2,:), bins1(iDay,2,:) ]= hist(bhv.stim_to_move{iDay,2}, 0:0.05:1.8);
    if ~any(h1(iDay,2,:))
        h1(iDay,2,:) = nan(1,37);
    end
    [h1(iDay,3,:), bins1(iDay,3,:) ]= hist(bhv.stim_to_move{iDay,3}, 0:0.05:1.8);
    if ~any(h1(iDay,3,:))
        h1(iDay,3,:) = nan(1,37);
    end
    end
    if any(h1(end,3,:)) 
        
         im = imagesc([squeeze(h1(:,1,:)), squeeze(h1(:,2,:)),squeeze(h1(:,3,:))])
         set(im, 'AlphaData', ~isnan(get(im, 'CData')));
         set(gca, 'color', [0.5, 0.5, 0.5]);
         xticks([1:18.5:37*3+1])
         xticklabels(strtrim(cellstr(num2str([0:0.9:1.8*3]'))'))
         
    else
        imagesc(squeeze(bins1(iDay,1,:))', [],squeeze(h1(:,1,:)))
    end

    colormap(brewermap([],'*RdBu'));
    caxis([- max(abs(max(max(im.CData)))),  max(abs(max(max(im.CData))))])
    catch
    end
    subplot(235)
%     plot(rxn_measured_med(:,1)); hold on;
%     plot(rxn_alt_med(:,1))
%     rxn_measured_allcat=[bhv.stim_to_move{:,1}];
%     
%     rxn_measured_med = accumarray([day_idx_go1'], ...
%     rxn_measured_allcat.*AP_nanout(rxn_measured_allcat < 0.1), ...
%     [size(h1,1)],@(x) nanmedian(x),NaN);
%     % nanmedian for each day ? 
    % permute 
%     n_rxn_altsample = 1000;
% rxn_alt = cellfun(@(rxn,use_days,use_trials) ...
%     cellfun(@(rxn,use_trials) ...
%     cell2mat(cellfun(@(x) datasample(x,n_rxn_altsample)',rxn(use_trials),'uni',false)), ...
%     rxn(use_days),use_trials(use_days),'uni',false), ...
%     {bhv.alt_stim_move_t},use_days,use_rxn,'uni',false)';

%     rxn_alt_med = cell2mat(permute(arrayfun(@(x) ...
%     accumarray([day_idx,animal_idx], ...
%     rxn_alt_allcat(:,x).*AP_nanout(rxn_alt_allcat(:,x) < 0.1), ...
%     [max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
%     1:n_rxn_altsample,'uni',false),[1,3,2]));
% % 
%     rxn_measured_med = accumarray([day_idx,animal_idx], ...
%     rxn_measured_allcat.*AP_nanout(rxn_measured_allcat < 0.1), ...
%     [max_days,length(bhv)],@(x) nanmedian(x),NaN);
% rxn_alt_med = cell2mat(permute(arrayfun(@(x) ...
%     accumarray([day_idx,animal_idx], ...
%     rxn_alt_allcat(:,x).*AP_nanout(rxn_alt_allcat(:,x) < 0.1), ...
%     [max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
%     1:n_rxn_altsample,'uni',false),[1,3,2]));


    transparencyValues = 0:1 / length(noGoDay):1;
    bhvOut(iAnimal).protocol_idx = protocol_idx;
    bhvOut(iAnimal).day_labels = day_labels;
    bhvOut(iAnimal).binBorders = binBorders;
%    bhvOut(iAnimal).movingFracGo1 = bhv.movingFrac;
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
    bhvOut(iAnimal).goLeft1 =  bhv.goLeft1;
    bhvOut(iAnimal).goLeft2 =  bhv.goLeft2;
    bhvOut(iAnimal).noGo1 =  bhv.noGo1;
    bhvOut(iAnimal).noGo2 =  bhv.noGo2;
    bhvOut(iAnimal).nTrials1 =  bhv.nTrials1;
    bhvOut(iAnimal).nTrials2 =  bhv.nTrials2;
    bhvOut(iAnimal).repeatOnMisses = bhv.repeatOnMisses;
    bhvOut(iAnimal).response_trials = bhv.response_trials;
    bhvOut(iAnimal).hitValues = bhv.hitValues;
    bhvOut(iAnimal).go1trials = bhv.go1trials;
    bhvOut(iAnimal).noGotrials = bhv.noGotrials;
    bhvOut(iAnimal).stim_to_move = bhv.stim_to_move;


    
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
            meanGo1Correct(iNoGoday, :) = nan(1600,1);
            catch
                meanGo1Correct(iNoGoday, :) = NaN;
            end
        meanGo2Correct(iNoGoday, :) = NaN;
        meanNoGoCorrect(iNoGoday, :) = nan(1600,1);
        meanGo1Incorrect(iNoGoday, :) = nan(1600,1);
        meanGo2Incorrect(iNoGoday, :) = nan(1600,1);
        meanNoGoIncorrect(iNoGoday, :) = nan(1600,1);
            

        end
    end

    an(iAnimal).noGo = bhv.noGo;
    an(iAnimal).nTrials = bhv.nTrials;
    an(iAnimal).goLeft = bhv.goLeft;
    an(iAnimal).noGoDay = noGoDay;

    an(iAnimal).stim_rxn_time = bhv.stim_rxn_time;
%    an(iAnimal).movingFrac = bhv.movingFrac;
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

            
            scatter([1, 2, 3], nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2])), [], col1, 'filled'); hold on;
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

clearvars meanGoLeftAcrossAnimals
clearvars meanNoGoAcrossAnimals
figure();
for iAnimal = 1:length(animalsAll)
    try
    transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
    %for iDay = 1:size(an(iAnimal).noGoDay, 2)
    col1 = rgb('DarkBlue');
    col2 = rgb('DarkRed');
    if ismember(iAnimal, 1:4)
        meanGoLeftAcrossAnimals(iAnimal,:) = nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
        meanNoGoAcrossAnimals(iAnimal,:) = nanmean(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
    else
        meanGoLeftAcrossAnimals(iAnimal,:) = 1 - nanmean(an(iAnimal).goLeft(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
        meanNoGoAcrossAnimals(iAnimal,:) = nanmean(an(iAnimal).noGo(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]))./nanmean(an(iAnimal).nTrials(an(iAnimal).noGoDay(end-nDays:end), [1, 3, 2]));
   
    end
    catch
    end
end
try
scatter([1, 2, 3], nanmean(meanGoLeftAcrossAnimals), [], col1, 'filled'); hold on;
p = plot([1, 2, 3], nanmean(meanGoLeftAcrossAnimals), 'Color', col1);
hold on;
errorbar([1, 2, 3],nanmean(meanGoLeftAcrossAnimals),...
    nanstd(meanGoLeftAcrossAnimals)./length(animalsAll), 'Color',col1)
makepretty;

scatter([1, 2, 3], nanmean(meanNoGoAcrossAnimals), [], col2, 'filled')
p = plot([1, 2, 3], nanmean(meanNoGoAcrossAnimals), 'Color', col2);

errorbar([1, 2, 3],nanmean(meanNoGoAcrossAnimals),...
    nanstd(meanNoGoAcrossAnimals)./length(animalsAll), 'Color',col2)
makepretty;
axis square

xlabel('image type')
ylabel('frac \color[rgb]{0,0,1}go left \color[rgb]{1,0,0} no go')
xticks([1, 2, 3])
xticklabels({'Go1', 'NoGo', 'Go2'})
ylim([0, 1])
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['mean +/- s.e. across ' num2str(length(animalsAll)) ' mice'], 'Box','off')
makepretty;
grid on;
catch
end


if plotAll
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


    if plotAll == true
        figure();

        for iAnimal = 5:6
            subplot(2, 1, iAnimal-4)
            transparencyValues = 0:1 / length(an(iAnimal).noGoDay):1;
            if iAnimal == 6
                thiD = size(an(iAnimal).noGoDay, 2) - 2:size(an(iAnimal).noGoDay, 2);
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
        if plotAll
            figure(); %last two days
            iAnimal = 6
            plot(nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), 'Color', rgb('RoyalBlue'))
            makepretty;
            hold on;
            errorbar([1, 2], nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), ...
                nanstd([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]./sqrt(2)), ...
                'Color', rgb('RoyalBlue'), 'LineWidth', 2)
            makepretty;
            hold on;
            plot(nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), 'Color', rgb('IndianRed'))
            makepretty;
            hold on;
            errorbar([1, 2], nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), ...
                nanstd([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]./sqrt(2)), ...
                'Color', rgb('IndianRed'), 'LineWidth', 2)
            makepretty;
            hold on;

            iAnimal = 5
            plot(nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), 'Color', rgb('Blue'))
            makepretty;
            hold on;
            errorbar([1, 2], nanmean([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), ...
                nanstd([an(iAnimal).goLeft(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]./sqrt(2)), ...
                'Color', rgb('Blue'), 'LineWidth', 2)
            makepretty;
            hold on;

            hold on;
            plot(nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), 'Color', rgb('Red'))
            makepretty;
            hold on;
            errorbar([1, 2], nanmean([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]), ...
                nanstd([an(iAnimal).noGo(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]) ./ (an(iAnimal).nTrials(animalsPhase4{iAnimal-4}(end, -1:end), [1, 3]))]./sqrt(2)), ...
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


        for iMouse = [5:6]
            subplot(2, 1, iMouse-4)
            transparencyValues = 0:1 / length(an(iMouse).noGoDay):1;
            for iDay = 1:size(an(iMouse).noGoDay, 2)
                if iDay == size(an(iMouse).noGoDay, 2)
                    col1 = rgb('DarkBlue');
                    col2 = rgb('DarkRed');

                else
                    col1 = 'b';
                    col2 = 'r';
                end
                scatter([1, 2, 3], an(iMouse).goLeft(an(iMouse).noGoDay(iDay), [1, 3, 2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1, 3, 2]), [], col1, 'filled', ...
                    'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
                p = plot([1, 2, 3], an(iMouse).goLeft(an(iMouse).noGoDay(iDay), [1, 3, 2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1, 3, 2]), 'Color', col1);
                p.Color(4) = transparencyValues(iDay+1);
                hold on;
                makepretty;
                scatter([1, 2, 3], an(iMouse).noGo(an(iMouse).noGoDay(iDay), [1, 3, 2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1, 3, 2]), [], col2, 'filled', ...
                    'MarkerFaceAlpha', transparencyValues(iDay+1), 'MarkerEdgeAlpha', transparencyValues(iDay+1))
                p = plot([1, 2, 3], an(iMouse).noGo(an(iMouse).noGoDay(iDay), [1, 3, 2])./an(iMouse).nTrials(an(iMouse).noGoDay(iDay), [1, 3, 2]), 'Color', col2);
                p.Color(4) = transparencyValues(iDay+1);
                makepretty;
            end

            xlabel('image type')
            ylabel('frac go left(blue)/nogo (red)')
            xticks([1, 2, 3])
            xticklabels({'Go1', 'NoGo', 'Go2'})
            ylim([0, 1])
            makepretty;
            grid on;
        end
    end
end
end