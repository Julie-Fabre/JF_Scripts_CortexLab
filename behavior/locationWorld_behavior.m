clear all;
clc;

% 2021 04 08 added no go analysis 
animalsAll = {'JF036', 'JF037', 'JF038', 'JF039'}; %JF029','JF032', 'JF033', 'JF034', 'JF035' };%{'JF032'};%{'JF025', 'JF026', 'JF028', 'JF029','JF032', 'JF033', 'JF034', 'JF035' };%chnage to the names of the animals you want here
for iAnimal = 1:size(animalsAll, 2)

    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'location'; %protocol name contains this name
    protocol2 = 'goNoGo';
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, protocol2); %find the experiments with the choiceworld protocol
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day = [];
    for curr_day = 1:length(experiments)

        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        trNum = [];
        % If multiple experiments, only use the largest one (usually multiple
        % happens if mess ups first/last one is good)
        for curr_experiment = 1:length(experiment_num)
            experiment = experiment_num(curr_experiment);

            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
            load(block_filename)

            trNum(curr_experiment) = length(block.events.newTrialTimes);
        end
        
        for curr_experiment = find(trNum == max(trNum))

            experiment = experiment_num(curr_experiment);

            [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');
            load(block_filename)
            %correct non repeat trials
            response_trials = 1:length(block.events.responseValues);
            repeatOnMisses = block.events.repeatTrialValues(response_trials) > 0;

            correctTrialsLeft = block.events.trialSideValues(response_trials) == -1 & block.events.hitValues(response_trials) == 1;
            correctTrialsRight = block.events.trialSideValues(response_trials) == 1 & block.events.hitValues(response_trials) == 1;
            sum(correctTrialsLeft);
            sum(correctTrialsRight);
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);


            response_trials = 1:length(block.events.responseValues);
            clearvars correctTrials
            iImg = 1;
            leftGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == -1));
            rightGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == 1));
            leftT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == -1));
            rightT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == 1));

            correctTrials((iImg - 1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & ...
                block.events.hitValues(response_trials) == 1);

            correctTrials((iImg - 1)*2+2) = sum(block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 ...
                );


            %                     bhv.n_trials_condition(curr_day, :) = performance(2, :);
            %                 bhv.go_left_trials(curr_day, :) = performance(end, :);
            %
            %                     correctTrials((iImg-1)*2+1) = sum(~repeatOnMisses(response_trials) & block.events.trialSideValues(response_trials) == -1 & block.events.hitValues(response_trials) == 1 ...
            %                         & block.events.stimNValues(response_trials) == imgs(iImg));
            %                      nTrials((iImg-1)*2+1) = sum(diff(rightT([response_trials,response_trials(end)+1])) ==1 &...
            %                          block.events.stimNValues(response_trials) == imgs(iImg));
            % %                       correctTrials((iImg-1)*2+2) = sum(~repeatOnMisses(response_trials) & block.events.trialSideValues(response_trials) == 1 & block.events.hitValues(response_trials) == 1 ...
            % %                         & block.events.stimNValues(response_trials) == imgs(iImg));
            %
            %                     nTrials((iImg-1)*2+2) = sum(diff(leftT([response_trials,response_trials(end)+1])) ==1 &...
            %                          block.events.stimNValues(response_trials) == imgs(iImg));

            nTrials((iImg - 1)*2+1) = sum(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == -1 ...
                );

            nTrials((iImg - 1)*2+2) = sum(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == 1 ...
                );

            conditi = [-1, 1];


            goLeft((iImg - 1)*2+1) = sum(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == -1 & ...
                block.events.hitValues(response_trials) == 0 ...
                );

            goLeft((iImg - 1)*2+2) = sum(~[repeatOnMisses(response_trials(2:end)), 0] & ...
                block.events.trialSideValues(response_trials) == 1 & ...
                block.events.hitValues(response_trials) == 1 ...
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
            win = [-5, 20];
            binSize=0.1;
            binBorders = win(1):binSize:win(2);
            binArray=[];
            for r = 1:length(block.events.stimOnTimes)-1
                  [n,binCenters] = histdiff(movingRN_times, block.events.stimOnTimes(r), binBorders);
                    binArray(r,:) = n;
                    %find(movingRN_times > block.events.stimOnTimes(r)-win(1) & movingRN_times < block.events.stimOnTimes(r+1))
            end
            cp=cumsum(binArray>1);
            movingFrac = cp(end, :) / size(binArray,1);

           
            %                     goLeft((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & ...
            %                         diff(leftGood([response_trials,response_trials(end)+1])) == 0 &...
            %                         block.events.stimNValues(response_trials) == imgs(iImg));
            %                     goLeft((iImg-1)*2+2) = sum(block.events.trialSideValues(response_trials) == 1 &...
            %                         diff(rightGood([response_trials,response_trials(end)+1])) > 0 &...
            %                         block.events.stimNValues(response_trials) == imgs(iImg));

            bhv.correctTrials(curr_day, :) = correctTrials;
            bhv.nTrials(curr_day, :) = nTrials;
            bhv.n_trials(curr_day, :) = n_trials;
            bhv.total_water(curr_day, :) = total_water;
            bhv.conditions(curr_day, :) = conditi;
            bhv.goLeft(curr_day, :) = goLeft;
            keep_day = [keep_day, curr_day];


            %stims
        end
    end

    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels = cellfun(@(day) [' ', day(6:end)], ...
        {experiments(keep_day).day}, 'uni', false);

    figure('Name', animal);
    subplot(211)
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
    set(gca, 'XTick', day_num);
    set(gca, 'XTickLabel', day_labels);
    set(gca, 'XTickLabelRotation', 90);
    makepretty;

    subplot(212)
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
    xlabel('< 0 = left, > 0 = right');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(keep_day));
    set(gca, 'YTickLabel', day_labels);

    for iDay = 1:size(keep_day, 2)
        txt = num2str(bhv.correctTrials(keep_day(iDay), :)');
        text((1:length(con))-0.2, ones(1, length(con))*iDay, txt, 'BackgroundColor', 'w')
    end
    %axis square;
    caxis([0, 1])
    makepretty;
    
    an(iAnimal).movingFrac = movingFrac; 

end

figure();
plot(binBorders(1:end-1), an(1).movingFrac);
hold on; 
makepretty;
plot(binBorders(1:end-1), an(2).movingFrac);
hold on; 
makepretty;
plot(binBorders(1:end-1), an(3).movingFrac);
hold on; 
makepretty;
plot(binBorders(1:end-1), an(4).movingFrac);
hold on; 
legend(animalsAll)
xlabel('time from stim onset (s)')
ylabel('fraction moving')
makepretty;