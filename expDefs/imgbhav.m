animalsAll = {'JF025', 'JF026'}; %JF029','JF032', 'JF033', 'JF034', 'JF035' };%{'JF032'};%{'JF025', 'JF026', 'JF028', 'JF029','JF032', 'JF033', 'JF034', 'JF035' };%chnage to the names of the animals you want here
for iAnimal = 1:size(animalsAll, 2)

    animal = animalsAll{1, iAnimal}; %this animal
    protocol = 'Imgs'; %protocol name contains this name
    protocol2 = 'Imgs';
    flexible_name = true; %protocol name can be slightly different
    experiments = AP_find_experimentsJF(animal, protocol, protocol2); %find the experiments with the choiceworld protocol
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day=[];
    for curr_day = 1:length(experiments)

        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;

        % If multiple experiments, only use the last one (usually multiple
        % happens if mess ups and final one is good)
        for curr_experiment = 1%length(experiment_num)

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
                    correctTrials((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & block.events.hitValues(response_trials) == 1 ...
                        & block.events.stimNValues(response_trials) == imgs(iImg));
                    nTrials((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));
                    correctTrials((iImg-1)*2+2) = sum(block.events.trialSideValues(response_trials) == 1 & block.events.hitValues(response_trials) == 1 ...
                        & block.events.stimNValues(response_trials) == imgs(iImg));
                    nTrials((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == 1 & ...
                        block.events.stimNValues(response_trials) == imgs(iImg));
                    conditi((iImg-1)*2+1) = -imgs(iImg);
                    conditi((iImg-1)*2+2) = imgs(iImg);
                    goLeft((iImg-1)*2+1) = sum(block.events.trialSideValues(response_trials) == -1 & block.events.hitValues(response_trials) == 0 ...
                        & block.events.stimNValues(response_trials) == imgs(iImg));
                    goLeft((iImg-1)*2+2) = sum(block.events.trialSideValues(response_trials) == 1 & block.events.hitValues(response_trials) == 1 ...
                        & block.events.stimNValues(response_trials) == imgs(iImg));
                end
                bhv.correctTrials(curr_day, :,:) = correctTrials; 
                bhv.nTrials(curr_day, :,:) = nTrials;
                bhv.n_trials(curr_day, :) = n_trials;
                bhv.total_water(curr_day, :) = total_water; 
                bhv.conditions(curr_day,:) = conditi;
                bhv.goLeft(curr_day,:) = goLeft;
                keep_day = [keep_day, curr_day];
            else
                
            end
            %stims
        end
    end
    
    day_num = cellfun(@(x) datenum(x), {experiments(keep_day).day});
    day_labels = cellfun(@(day) [ ' ', day(6:end)], ...
        {experiments(keep_day).day}, 'uni', false);
    
    figure('Name', animal);
    subplot(211)
    yyaxis left
    
    scatter(day_num, bhv.n_trials(keep_day,:)); hold on;
    plot(day_num, bhv.n_trials(keep_day,:), 'linewidth', 2);
    ylabel('Trials');
    yyaxis right
    
    scatter(day_num, bhv.total_water(keep_day,:), 'linewidth', 2); hold on;
    plot(day_num, bhv.total_water(keep_day,:), 'linewidth', 2);
    ax = gca;
    hold on;

    ylabel('Total water (ul)');
    xlabel('Session');
    set(gca, 'XTick', day_num);
    set(gca, 'XTickLabel', day_labels);
    set(gca, 'XTickLabelRotation', 90);
    makepretty;
    
    subplot(212)
    con = bhv.conditions(keep_day,:);
    ccc = num2str(con);
    cc  = textscan(ccc,'%s','Delimiter',' ');
    cc=cc{1};
    cc(strcmp('',cc)) = [];
    im = imagesc(1:length(con), 1:size(bhv.correctTrials(keep_day,:)), bhv.goLeft(keep_day,:)./bhv.nTrials(keep_day,:));
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));
    set(gca, 'color', [0.5, 0.5, 0.5]);
    colormap(brewermap([], '*RdBu'));
    c = colorbar;
    ylabel(c, 'Go left (frac)');
    xlabel('Condition');
    ylabel('Session');
    set(gca, 'XTick', 1:length(con));
    set(gca, 'XTickLabel', cc);
    set(gca, 'YTick', 1:length(experiments));
    set(gca, 'YTickLabel', day_labels);
    txt = num2str(squeeze(bhv.correctTrials(keep_day,1,:)));
    text((1:length(con))-0.2,ones(1,length(con)),txt, 'BackgroundColor', 'w')
    %axis square;
    caxis([0, 1])
    makepretty; 
    
    figure();
    for iImg = 1:8
        subplot(1,8,iImg )
        thisIm = imread(['\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img' num2str(iImg), '.jpeg']);
        imagesc(thisIm)
        colormap(gray)
        axis square; 
        title(num2str(iImg))
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'YTick',[]);
        set(gca,'XTick',[]);
        makepretty; 
    end
    
end