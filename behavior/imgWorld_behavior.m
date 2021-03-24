clear all;
clc;

animalsAll = {'JF026', 'JF028'}; %JF029','JF032', 'JF033', 'JF034', 'JF035' };%{'JF032'};%{'JF025', 'JF026', 'JF028', 'JF029','JF032', 'JF033', 'JF034', 'JF035' };%chnage to the names of the animals you want here
for iAnimal = 1:size(animalsAll, 2)

    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'Imgs'; %protocol name contains this name
    protocol2 = 'Imgs';
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, protocol2); %find the experiments with the choiceworld protocol
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day = [];
    for curr_day = length(experiments)

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
        % If multiple experiments, only use the last one (usually multiple
        % happens if mess ups and final one is good)
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

            if isfield(block.events, 'stimNValues')
                imgs = unique(block.events.stimNValues);
                response_trials = 1:length(block.events.stimNValues);
                clearvars correctTrials
                for iImg = 1:length(imgs)
                    hit = [block.events.hitValues(2:end), 0];
            	    miss = [block.events.missValues(2:end), 0];
                    leftGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == -1));
                    rightGood = block.events.sessionPerformanceValues(3, find(block.events.sessionPerformanceValues(1, :) == 1));
                    leftT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == -1));
                    rightT = block.events.sessionPerformanceValues(2, find(block.events.sessionPerformanceValues(1, :) == 1));

                    correctTrials((iImg - 1)*2+1) = sum( block.events.trialSideValues(response_trials) == -1 & ...
                        hit(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));

                    correctTrials((iImg - 1)*2+2) = sum( block.events.trialSideValues(response_trials) == 1 & ...
                        hit(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));


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

                    nTrials((iImg - 1)*2+1) = sum( ...
                        block.events.trialSideValues(response_trials) == -1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));

                    nTrials((iImg - 1)*2+2) = sum( ...
                        block.events.trialSideValues(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));
                    truC = [~[repeatOnMisses(response_trials(2:end)), 0] & block.events.trialSideValues(response_trials) == 1 & ...
                        hit(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg)]*2 - [~[repeatOnMisses(response_trials(2:end)), 0] & ...
                        block.events.trialSideValues(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg)];

                    conditi((iImg - 1)*2+1) = -imgs(iImg);
                    conditi((iImg - 1)*2+2) = imgs(iImg);

                    goLeft((iImg - 1)*2+1) = sum( ...
                        block.events.trialSideValues(response_trials) == -1 & ...
                        miss(response_trials) == 0 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));

                    goLeft((iImg - 1)*2+2) = sum( ...
                        block.events.trialSideValues(response_trials) == 1 & ...
                        hit(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));

                    %                     goLeft((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & ...
                    %                         diff(leftGood([response_trials,response_trials(end)+1])) == 0 &...
                    %                         block.events.stimNValues(response_trials) == imgs(iImg));
                    %                     goLeft((iImg-1)*2+2) = sum(block.events.trialSideValues(response_trials) == 1 &...
                    %                         diff(rightGood([response_trials,response_trials(end)+1])) > 0 &...
                    %                         block.events.stimNValues(response_trials) == imgs(iImg));
                end
                bhv.correctTrials(curr_day, :) = correctTrials;
                bhv.nTrials(curr_day, :) = nTrials;
                bhv.n_trials(curr_day, :) = n_trials;
                bhv.total_water(curr_day, :) = total_water;
                bhv.conditions(curr_day, :) = conditi;
                bhv.goLeft(curr_day, :) = goLeft;
                keep_day = [keep_day, curr_day];
            else

            end
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
    con = bhv.conditions(keep_day, :);
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
    xlabel('Image #, < 0 = left, > 0 = right');
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

    an(iAnimal).bhv = bhv;
    an(iAnimal).keep_day = keep_day;
end

figure();
for iImg = 1:8
    subplot(2, 8, iImg)
    thisIm = imread(['\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7StimsNew\img', num2str(iImg), '.jpeg']);
    imagesc(thisIm)
    colormap(gray)
    axis square;
    title(num2str(iImg))
    set(gca, 'YTickLabel', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    makepretty;
end
subplot(2, 8, [9:16])
%percentage correct
colorsT = {'blue','bo';'black','ko'; 'green', 'go'};
for iAnimal = 1:size(an, 2)
    correctTrialsPerStim = reshape(an(iAnimal).bhv.correctTrials(), size(an(iAnimal).bhv.correctTrials, 1)*2, size(an(iAnimal).bhv.correctTrials, 2)/2);
    nTrialsPerStim = reshape(an(iAnimal).bhv.nTrials, size(an(iAnimal).bhv.nTrials, 1)*2, size(an(iAnimal).bhv.nTrials, 2)/2);
    val(iAnimal,:) = sum(correctTrialsPerStim, 1)./sum(nTrialsPerStim, 1);
    scatter(1:8, val(iAnimal,:), 10, colorsT{iAnimal,1}, 'filled');
    hold on;
    makepretty;
end

errorbar(1:8, nanmean(val),nanstd(val),'ro','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2)
    hold on;
    makepretty;
    
    for iAnimal = 1:size(an, 2)
    correctTrialsPerStim = reshape(an(iAnimal).bhv.correctTrials(), size(an(iAnimal).bhv.correctTrials, 1)*2, size(an(iAnimal).bhv.correctTrials, 2)/2);
    nTrialsPerStim = reshape(an(iAnimal).bhv.nTrials, size(an(iAnimal).bhv.nTrials, 1)*2, size(an(iAnimal).bhv.nTrials, 2)/2);
    val(iAnimal,:) = sum(correctTrialsPerStim, 1)./sum(nTrialsPerStim, 1);
    scatter(1:8, val(iAnimal,:), 30, colorsT{iAnimal,1}, 'filled');
    hold on;
    makepretty;
    end
xlabel('image #')
ylabel('percentage correct trials')
    legend({animalsAll{:},'mean +/- std'})
    makepretty;