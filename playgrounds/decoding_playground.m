%%

% to run glmnet, had to :
% - run `mex -setup fortran` - % actually was unnecessary - the compiled library was alread present 
% - install libgfortran3 :  sudo apt-get install libgfortran3. I had gcc
% version 8.4.0 installed when I did this. 

% glmnet() logistic regression
% can use either l1 penalty or just explicitly select e.g. top % of cells
% (highest difference score) ]
% statistics experiment by experiment 
% question 1: DMS, easier to decode stim 3 than 1+2 in go no go?
% question 2: DMS, hard to decode 1+2+3 in go go go?
keep task_data_here
passive = 1;
if passive
task_data_here = load('C:\Users\julie\Dropbox\MATLAB\task_data_passive.mat');
idx = 5; %2: task, 4/5 : passive

else
task_data_here = load('C:\Users\julie\Dropbox\MATLAB\task_data_goNogo3.mat');
idx = 2; %2: task, 4/5 : passive

end

plot_regions = [1,2,3];%[1, 2, 5]; % Striatum, GPe, SNr 
for iSession = 1:size(task_data_here.av_per_trial,2)
    session_nCells(iSession) = size(task_data_here.av_per_trial{idx,iSession},1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells>0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);

%% fraction of cells, average response
% per session

for iRegion = 1:size(plot_regions,2) 
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);% & ...
        %task_data_here.pvalue_shuffled_005{1,idx}' == 1;
    these_units_pval = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2) & ...
        task_data_here.pvalue_shuffled_005{1,idx}' == 1;
   unitCount = 0;
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        trials_no_move = task_data_here.trial_types{idx,thisSession}(task_data_here.no_move_trials{idx,thisSession},:);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));
        
        
        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move, trialTypes, 'rows');%QQ
            theseTrialTypes = trials_no_move(theseTrials,:);
            these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
            these_units_session_pval = these_units_pval(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);

            activity_per_trial_neuron = task_data_here.av_per_trial{idx,thisSession}(these_units_session,theseTrials)';

            n_cells_session(iRegion,iSession) = sum(these_units_session);
            %n_cells_pval_session(iRegion,iSession) = sum(these_units_session_pval);

               activity_per_trial_neuron = task_data_here.av_per_trial{idx,thisSession}(these_units_session,theseTrials)';
              for_baseline_per_neuron_per_condition = task_data_here.av_psth{idx,thisSession}(these_units_session,:,1:200);


               if size(for_baseline_per_neuron_per_condition,2) == 39
                  cond_inds = [10,34; 16,34; 16, 10];
               elseif size(for_baseline_per_neuron_per_condition,2) == 26
                  cond_inds = [7,23; 11,23; 11, 7];
               elseif size(for_baseline_per_neuron_per_condition,2) == 4
                  cond_inds = [2,4; 3,4; 3,2];
               else
                  disp('wtf')
               end
               clearvars is_pval
    ff= find(these_units_session);
               if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron >0)))
                 for iNeuron = 1:size(activity_per_trial_neuron,2)
                   for iPair = 1:3 
                       theseTrials_pval = ismember(trials_no_move, trialTypes(iPair,:), 'rows');%QQ
                       [p,h] = signrank(task_data_here.av_per_trial_base{idx,thisSession}(ff(iNeuron),theseTrials_pval),...
                           task_data_here.av_per_trial{idx,thisSession}(ff(iNeuron),theseTrials_pval));
                       is_pval(iPair,iNeuron) = p <= 0.05;
                   end
                 end
                 for iPair = 1:3
                     n_cells_pval_session(iRegion,iPair,iSession) = sum(is_pval(iPair,:));
                 end
               end
               
        end
    end
end

figure(); hold on;
cols = lines(3);
for iRegion = 1:3
    subplot(3,1,iRegion); hold on;
    for iPair =1:3
    pct_per_session = [squeeze(n_cells_pval_session(iRegion,iPair,:))' ./ n_cells_session(iRegion,:) .*100];
    scatter(ones(size(pct_per_session,2),1).*iPair,pct_per_session, 15,cols(iPair,:),'filled')
    bar(iPair, nanmean(pct_per_session ),'FaceColor',cols(iPair,:))
    end
end
ylabel('% visual cells')

figure();


for iRegion = 1:size(plot_regions,3)
    these_units = task_data_here.unit_area == iRegion;% & (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);
   
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession);
        trials_no_move = task_data_here.trial_types{idx,thisSession}(task_data_here.no_move_trials{idx,thisSession},:);
        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move,...
                trialTypes, 'rows');%QQ
            theseTrialTypes = trials_no_move(theseTrials,:);
            these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1));

            
            if ~isempty(these_units_session) && ~islogical(these_units_session)
                figure();
            imagesc(squeeze(task_data_here.psth{1,2}(these_units_session,3,11,:)))
                activity_per_trial_neuron = task_data_here.av_per_trial{idx,thisSession}(:,theseTrials)';

                family = 'multinomial';%'binomial'; %logistic
               % fit = glmnet(activity_per_trial_neuron, theseTrialTypes, family, options);
               if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron >0)))
                   % get a good lambda value (?) 
                   
                    folds_trials = randsample(5, size(theseTrialTypes,1), true);

                    fit_cv_lambda = cvglmnet(activity_per_trial_neuron(:,:),...
                        theseTrialTypes(:,1), family, [], [],  5);

                    use_lambda = fit_cv_lambda.lambda_min;
                    fit_cv_lambda.glmnet_fit

                    %cvglmnetPlot(fit_cv_lambda);
                    %cvglmnetPredict(fit_cv_lambda,activity_per_trial_neuron(:,:),'lambda_min')
                    
                    
                    for iFold = 1:5 %divide into 5, train on 4/5th and test on remaining 1/5th

                        theseFolds = ones(5,1);
                        theseFolds(iFold) = 0;
                        options.lambda = use_lambda;
                        fit_cv = glmnet(activity_per_trial_neuron(ismember(folds_trials, find(theseFolds)),:), ...
                            theseTrialTypes(ismember(folds_trials, find(theseFolds)),1), family, options);
                        %glmnetPredict(object, newx, s, type, exact, offset)
                        prediction_cv = glmnetPredict(fit_cv, activity_per_trial_neuron(ismember(folds_trials, iFold),:), [], [], true);%, use_lambda, 'link');
                        [~ , trials_pred_thisType] = min(prediction_cv');
                        for iTrialType = 1:3
                            trials_thisType = double(theseTrialTypes(ismember(folds_trials, iFold),1) == trialTypes(iTrialType,1)).*iTrialType;
                            
                            accuracy_cv(iFold, iTrialType) = sum(trials_thisType==trials_pred_thisType')./sum(trials_thisType); %which trials correctly labelled
                        end

                    
                    end
                    disp(accuracy_cv)
                    figure();
                    subplot(121)
                    imagesc(accuracy_cv)
                    colorbar;
                    subplot(122)
                    
               end
            end
        end

    end
    %g2=randsample(2,100,true);
    %x=randn(100,20);
    %fit2=glmnet(x,g2,'binomial');
end

%% d' against shuffle. get sum d' for each session. compare against distribution -> p value. combine p values with Fisher method (P-values) or Stouffer (Z scores)
clearvars d_prime pooled_sd_all ci
thesePairs = [1,2;1,3;2,3];

unitCount = 0;

for iRegion = 1:size(plot_regions,2) 
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);% & ...
        %task_data_here.pvalue_shuffled_005{1,idx}' == 1;
   unitCount = 0;
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        trials_no_move = task_data_here.trial_types{idx,thisSession}(task_data_here.no_move_trials{idx,thisSession},:);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));
        
        
        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move, trialTypes, 'rows');%QQ
            theseTrialTypes = trials_no_move(theseTrials,:);
           these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
             %figure();
             %i=2;
             %imagesc(squeeze(task_data_here.av_psth{idx,thisSession}(:,i,:)))
               activity_per_trial_neuron = task_data_here.av_per_trial{idx,thisSession}(these_units_session,theseTrials)';
               for_baseline_per_neuron_per_condition = task_data_here.av_psth{idx,thisSession}(these_units_session,:,1:200);
               if size(for_baseline_per_neuron_per_condition,2) == 39
                  cond_inds = [10,34; 16,34; 16, 10];
               elseif size(for_baseline_per_neuron_per_condition,2) == 26
                  cond_inds = [7,23; 11,23; 11, 7];
               elseif size(for_baseline_per_neuron_per_condition,2) == 4
                  cond_inds = [2,4; 3,4; 3,2];
               else
                  disp('wtf')
               end
                all_trials = ismember(theseTrialTypes, trialTypes, 'rows');
               if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron >0)))
                 for iNeuron = 1:size(activity_per_trial_neuron,2)
                   for iPair = 1:3 
                       trials_1 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair,1),:), 'rows');
                       trials_2 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair,2),:), 'rows');
                       
                            baseline_1 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair,1),:)./0.001);
                      
                            baseline_2 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair,2),:)./0.001);
                       if baseline_1 > 0.5
                       average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ? 
                       average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001)); 

                       sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001);%./sqrt(sum(trials_1));
                       sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001);%./sqrt(sum(trials_2));

                       pooled_sd = sqrt((sd_stim_1*sd_stim_1 + sd_stim_2*sd_stim_2)./2) + 1;
                       pooled_sd_all{iRegion}(iNeuron + unitCount, iPair) = pooled_sd;
                       d_prime{iRegion}(iNeuron + unitCount, iPair) = abs( average_stim_1 - average_stim_2)/pooled_sd;
                      
                       end
                   end
                    if baseline_1 > 0.5
                    d_prime_session{iRegion}(iNeuron) =  d_prime{iRegion}(iNeuron + unitCount, 1) -  d_prime{iRegion}(iNeuron + unitCount,2);
                       for iShuffle = 1:1000
                           trials_1 = randsample(1:size(all_trials,1),round(size(all_trials,1)./3));
                           trials_2 = randsample(1:size(all_trials,1),round(size(all_trials,1)./3));
                           trials_3 = randsample(1:size(all_trials,1),round(size(all_trials,1)./3));
                           average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ? 
                       average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001));
                       average_stim_3 = (nanmean(activity_per_trial_neuron(trials_3, iNeuron)./0.001));

                       sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001);%./sqrt(sum(trials_1));
                       sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001);%./sqrt(sum(trials_2));
                       sd_stim_3 = nanstd(activity_per_trial_neuron(trials_3, iNeuron)./0.001);%./sqrt(sum(trials_2));

                       pooled_sd_1 = sqrt((sd_stim_1*sd_stim_1 + sd_stim_2*sd_stim_2)./2) + 1;
                       pooled_sd_2 = sqrt((sd_stim_1*sd_stim_1 + sd_stim_3*sd_stim_3)./2) + 1;
                       %pooled_sd_all{iRegion}(iNeuron + unitCount, iPair) = pooled_sd;
                       d_prime_shuffled_1= abs( average_stim_1 - average_stim_2)/pooled_sd_1;
                       d_prime_shuffled_2= abs( average_stim_1 - average_stim_3)/pooled_sd_2;
                        d_prime_shuffled(iShuffle,iNeuron) = d_prime_shuffled_1 - d_prime_shuffled_2;

                       end
                       ci{iRegion}(iNeuron + unitCount, iPair) = (average_stim_1 - average_stim_2)/(average_stim_1 + average_stim_2 + 0.1);

                       %d_prime_z(iNeuron + unitCount, iPair, iRegion) = zscore((activity_per_trial_neuron(trials_1, iNeuron)))...
                       %    - ztrans((activity_per_trial_neuron(trials_2, iNeuron)));
                    end
                   
                 end

                 % figure();
                 % for iPair=1:3
                 %    subplot(3,1,iPair)
                 %    histogram(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes 
                 %    nanmedian(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair))
                 %    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2))])
                 % end
                 unitCount = unitCount + size(activity_per_trial_neuron,2);

               end
           
        end

    end
    % prctile(sum(abs(d_prime_shuffled),2),1:100)
    session_pvalue_dprime(iRegion,iSession) = sum(abs(d_prime_session{iRegion})) >= prctile(sum(abs(d_prime_shuffled),2),[95:0.001:100]) ||...
         sum(abs(d_prime_session{iRegion})) <= prctile(sum(abs(d_prime_shuffled),2),5);

    %PERCENTRANK = @(YourArray, TheProbes) reshape( mean( bsxfun(@le, YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) )

    figure();
    histogram(sum(abs(d_prime_shuffled),2)); hold on;
    ylims = ylim;
    line([sum(abs(d_prime_session{iRegion})),sum(abs(d_prime_session{iRegion}))], [ylims(1), ylims(2)])
end
 



%% d ' histograms 
clearvars d_prime pooled_sd_all ci
thesePairs = [1,2;1,3;2,3];

unitCount = 0;

for iRegion = 1:size(plot_regions,2) 
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);% & ...
        %task_data_here.pvalue_shuffled_005{1,idx}' == 1;
   unitCount = 0;
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        trials_no_move = task_data_here.trial_types{idx,thisSession}(task_data_here.no_move_trials{idx,thisSession},:);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));
        
        
        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move, trialTypes, 'rows');%QQ
            theseTrialTypes = trials_no_move(theseTrials,:);
           these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
             %figure();
             %i=2;
             %imagesc(squeeze(task_data_here.av_psth{idx,thisSession}(:,i,:)))
               activity_per_trial_neuron = task_data_here.av_per_trial{idx,thisSession}(these_units_session,theseTrials)';
               for_baseline_per_neuron_per_condition = task_data_here.av_psth{idx,thisSession}(these_units_session,:,1:200);
               if size(for_baseline_per_neuron_per_condition,2) == 39
                  cond_inds = [10,34; 16,34; 16, 10];
               elseif size(for_baseline_per_neuron_per_condition,2) == 26
                  cond_inds = [7,23; 11,23; 11, 7];
               elseif size(for_baseline_per_neuron_per_condition,2) == 4
                  cond_inds = [2,4; 3,4; 3,2];
               else
                  disp('wtf')
               end

               if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron >0)))
                 for iNeuron = 1:size(activity_per_trial_neuron,2)
                   for iPair = 1:3 
                       trials_1 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair,1),:), 'rows');
                       trials_2 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair,2),:), 'rows');
                       
                            baseline_1 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair,1),:)./0.001);
                      
                            baseline_2 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair,2),:)./0.001);
                       if baseline_1 > 0.5
                       average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ? 
                       average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001)); 

                       sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001);%./sqrt(sum(trials_1));
                       sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001);%./sqrt(sum(trials_2));

                       pooled_sd = sqrt((sd_stim_1*sd_stim_1 + sd_stim_2*sd_stim_2)./2) + 1;
                       pooled_sd_all{iRegion}(iNeuron + unitCount, iPair) = pooled_sd;
                       d_prime{iRegion}(iNeuron + unitCount, iPair) = abs( average_stim_1 - average_stim_2)/pooled_sd;
                       ci{iRegion}(iNeuron + unitCount, iPair) = (average_stim_1 - average_stim_2)/(average_stim_1 + average_stim_2 + 0.1);

                       %d_prime_z(iNeuron + unitCount, iPair, iRegion) = zscore((activity_per_trial_neuron(trials_1, iNeuron)))...
                       %    - ztrans((activity_per_trial_neuron(trials_2, iNeuron)));
                       end
                   end
                 end

                 % figure();
                 % for iPair=1:3
                 %    subplot(3,1,iPair)
                 %    histogram(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes 
                 %    nanmedian(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair))
                 %    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2))])
                 % end
                 unitCount = unitCount + size(activity_per_trial_neuron,2);

               end
           
        end

    end

end
 


% plot d prime
figure();
for iRegion=1:size(plot_regions,2)
for iPair=1:3
    subplot(3,size(plot_regions,2),iPair + (iRegion-1)*(size(plot_regions,2)))
    % if 0, remove
   %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
    kp = find( abs(d_prime{iRegion}(:,iPair)) ~= 0 & ~isinf( abs(d_prime{iRegion}(:,iPair))));
   
    histogram(d_prime{iRegion}(kp,iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes 
    %nanmedian(d_prime(kp,iPair,iRegion))
    %nanmean(d_prime(kp,iPair,iRegion))
    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime{iRegion}(kp,iPair))))])
    xlabel('d-prime')
    ylabel('# of neurons')
    %xlim([-1.45, 1.45])
    %makepretty_lite;
end
end


% plot sel. idx 
figure();
for iRegion=1:size(plot_regions,2)
for iPair=1:3
    subplot(3,size(plot_regions,2),iPair + (iRegion-1)*(size(plot_regions,2)))
    % if 0, remove
   %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
    kp = find( abs(ci{iRegion}(:,iPair)) ~= 0 & ~isinf( abs(ci{iRegion}(:,iPair))));
   
    histogram(ci{iRegion}(kp,iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes 
    %nanmedian(d_prime(kp,iPair,iRegion))
    %nanmean(d_prime(kp,iPair,iRegion))
    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(ci{iRegion}(kp,iPair))))])
    xlabel('Sel. idx')
    ylabel('# of neurons')
    %xlim([-1.45, 1.45])
    %makepretty_lite;
end
end


figure();
for iRegion=1:size(plot_regions,2)
%for iPair=1:3
    subplot(size(plot_regions,2),1,iRegion)
    % if 0, remove
  % kp = find(abs(d_prime(:,iPair,iRegion)) < 1.4141 & abs(d_prime(:,iPair,iRegion)) ~= 0);
  %  kp = find( abs(d_prime(:,iPair,iRegion)) ~= 0);
   
    scatter(squeeze(d_prime{iRegion}(:,1)),squeeze(d_prime{iRegion}(:,2)),4,'filled'); hold on; % 1.4142   -1.4142 = one of stims has 0 spikes 
    %nanmedian(d_prime(kp,iPair,iRegion))
    %nanmean(d_prime(kp,iPair,iRegion))
    %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime(kp,iPair,iRegion))))])
    xlabel('1 vs 2')
    ylabel('1 vs 3')
    %xlim([-1.4, 1.4])
   % ylim([-1.4, 1.4])
   

    set(gca,'DataAspectRatio',[1 1 1])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
     xlimits = xlim;
    ylimits = ylim;
    maxLim = max([xlimits, ylimits]);
    minLim = min([xlimits, ylimits]);
    line([minLim, maxLim], [minLim, maxLim], 'Color', rgb('Red'));
    %makepretty;
    %get()
    
    
%end
end

figure();
for iRegion=1:size(plot_regions,2)
%for iPair=1:3
    subplot(size(plot_regions,2),1,iRegion)
    % if 0, remove
  % kp = find(abs(d_prime(:,iPair,iRegion)) < 1.4141 & abs(d_prime(:,iPair,iRegion)) ~= 0);
  %  kp = find( abs(d_prime(:,iPair,iRegion)) ~= 0);
   
    
scatterHistDiff(squeeze(d_prime{iRegion}(:,1)),squeeze(d_prime{iRegion}(:,2)), '','', rgb('Blue'), 0.5)
%scatter(squeeze(d_prime{iRegion}(:,1)),squeeze(d_prime{iRegion}(:,2)),4,'filled'); hold on; % 1.4142   -1.4142 = one of stims has 0 spikes 
    %nanmedian(d_prime(kp,iPair,iRegion))
    %nanmean(d_prime(kp,iPair,iRegion))
    %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime(kp,iPair,iRegion))))])
    xlabel('1 vs 2')
    ylabel('1 vs 3')
    %xlim([-1.4, 1.4])
   % ylim([-1.4, 1.4])
   

    set(gca,'DataAspectRatio',[1 1 1])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
     xlimits = xlim;
    ylimits = ylim;
    maxLim = max([xlimits, ylimits]);
    minLim = min([xlimits, ylimits]);
    line([minLim, maxLim], [minLim, maxLim], 'Color', rgb('Red'));
    %makepretty;
    %get()
    
    
%end
end

scatterHistDiff(x, y, xeb, yeb, colors, histogramBinSize)

%% dot products: 1 v 2, 1 v 3, 2 v 3
thesePairs = [1,2;1,3;2,3];
clearvars dot_prod_per_session
% dot(A,B) 
clearvars averages_flat
clearvars xcorr_per_session
averages_flat{3,3} =[];
for iRegion = 1:size(plot_regions,2)
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2) & ...
        task_data_here.pvalue_shuffled_005{1,idx}' == 1;
        unitCount= 0;
   
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
        average_per_neuron = task_data_here.av_psth{idx,thisSession}(these_units_session,:,:);
        
        if size(for_baseline_per_neuron_per_condition,2) == 39
            cond_inds = [10,34; 16,34; 16, 10];
        elseif size(for_baseline_per_neuron_per_condition,2) == 26
            cond_inds = [7,23; 11,23; 11, 7];
        elseif size(for_baseline_per_neuron_per_condition,2) == 4
                  cond_inds = [2,4; 3,4; 3,2];
        else
            disp('wtf')
        end
        if size(average_per_neuron,1) == 1
            averages_flat{iRegion,1} = [ averages_flat{iRegion,1}; ...
            squeeze(average_per_neuron(:,cond_inds(1,1),: ))'];
        averages_flat{iRegion,2} = [ averages_flat{iRegion,2}; ...
            squeeze(average_per_neuron(:,cond_inds(1,2),: ))'];
        averages_flat{iRegion,3} = [ averages_flat{iRegion,3}; ...
            squeeze(average_per_neuron(:,cond_inds(2,1),: ))'];

        else
            averages_flat{iRegion,1} = [ averages_flat{iRegion,1}; ...
            squeeze(average_per_neuron(:,cond_inds(1,1),: ))];
        averages_flat{iRegion,2} = [ averages_flat{iRegion,2}; ...
            squeeze(average_per_neuron(:,cond_inds(1,2),: ))];
        averages_flat{iRegion,3} = [ averages_flat{iRegion,3}; ...
            squeeze(average_per_neuron(:,cond_inds(2,1),: ))];
        end
        
        % if iRegion == 3 && ~isempty(these_units_session)
        % figure();
        % subplot(311)
        % imagesc(squeeze(average_per_neuron(:,cond_inds(1,1),: )))
        % subplot(312)
        % imagesc(squeeze(average_per_neuron(:,cond_inds(1,2),: )))
        % subplot(313)
        % imagesc(squeeze(average_per_neuron(:,cond_inds(2,1),:)))
        % end
        for iPair=1:3
            average_stim_1 = smoothdata(squeeze(average_per_neuron(:, cond_inds(iPair,1) ,:)), 2,'gaussian', [50 0]);
            average_stim_2 =  smoothdata(squeeze(average_per_neuron(:, cond_inds(iPair,2) ,:)), 2,'gaussian', [50 0]);
% figure();
%             plot(squeeze(average_per_neuron(1, cond_inds(iPair,1) ,:))); hold on;
%             plot(average_stim_1(1,:))
            dot_prod_per_session{iRegion}(:,iPair,iSession) = dot(average_stim_1, average_stim_2);
            for ineuron = 1:size(average_stim_1,2)
                cc = corrcoef(average_stim_1(:,ineuron), average_stim_2(:,ineuron));
                xcorr_per_session{iRegion}(ineuron,iPair,iSession) = cc(2);
            end
        end



    end
end

figure();
for iRegion=1:size(plot_regions,2)
    for iPair=1:3
        subplot(3,1,iRegion); hold on;%size(plot_regions,2),iPair + (iRegion-1)*(size(plot_regions,2)))
        plot(task_data_here.t,nanmean(dot_prod_per_session{iRegion}(:,iPair,:),3))
        %makepretty_lite;
        %ylim([0.65, 0.8])
        makepretty;
        xlim([-0.2, 0.6])
    end
end

figure();
for iRegion=1:size(plot_regions,2)
    for iPair=1:3
        subplot(3,1,iRegion); hold on;%size(plot_regions,2),iPair + (iRegion-1)*(size(plot_regions,2)))
        plot(task_data_here.t,nanmean(xcorr_per_session{iRegion}(:,iPair,:),3))
        %makepretty_lite;
        %ylim([0.65, 0.8])
        makepretty;
        xlim([-0.15, 0.6])
    end
end

figure();
for iRegion=1:size(plot_regions,2)
    for iPair=1:3
        subplot(3,3,(iRegion-1)*3+iPair);
        imagesc(task_data_here.t,[],(averages_flat{iRegion,iPair}-nanmean(averages_flat{iRegion,iPair}(:,1:200),2))...
            ./(nanstd(averages_flat{iRegion,iPair}(:,1:200),[],2)))

    end
end


%% correlation (each neuron)
thesePairs = [1,2;1,3;2,3];
clearvars dot_prod_per_session
% dot(A,B) 
clearvars averages_flat
clearvars xcorr_per_session
averages_flat{3,3} =[];
xcorr_per_neuron_session{1,3} =[];
xcorr_per_neuron_session{2,3} =[];
xcorr_per_neuron_session{3,3} =[];
xcorr_per_neuron_session{1,2} =[];
xcorr_per_neuron_session{2,2} =[];
xcorr_per_neuron_session{3,2} =[];
xcorr_per_neuron_session{1,1} =[];
xcorr_per_neuron_session{2,1} =[];
xcorr_per_neuron_session{3,1} =[];
for iRegion = 1:size(plot_regions,2)
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);% & ...
        %task_data_here.pvalue_shuffled_005{1,2}' == 1;
        unitCount= 0;
   
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
        average_per_neuron = task_data_here.av_psth{2,thisSession}(these_units_session,:,:);
        
        if size(for_baseline_per_neuron_per_condition,2) == 39
            cond_inds = [10,34; 16,34; 16, 10];
        elseif size(for_baseline_per_neuron_per_condition,2) == 26
            cond_inds = [7,23; 11,23; 11, 7];
        else
            disp('wtf')
        end

       
        for iPair=1:3
            average_stim_1 = squeeze(average_per_neuron(:, cond_inds(iPair,1) ,:));
            average_stim_2 = squeeze(average_per_neuron(:, cond_inds(iPair,2) ,:));
            dot_prod_per_session{iRegion}(:,iPair,iSession) = dot(average_stim_1, average_stim_2);
            for ineuron = 1:size(average_stim_1,1)
                if size(average_stim_1,2) == 1
                     cc = corrcoef(average_stim_1(:), average_stim_2(:));
                xcorr_per_neuron_session{iRegion,iPair} = [xcorr_per_neuron_session{iRegion,iPair}, cc(2)];
                else
                     cc = corrcoef(average_stim_1(ineuron,:), average_stim_2(ineuron,:));
                    xcorr_per_neuron_session{iRegion,iPair} = [xcorr_per_neuron_session{iRegion,iPair}, cc(2)];
                end
               
            end
        end



    end
end



figure();
for iRegion=1:size(plot_regions,2)
    %for iPair=1:3
        subplot(3,1,iRegion); hold on;%size(plot_regions,2),iPair + (iRegion-1)*(size(plot_regions,2)))
        scatter(xcorr_per_neuron_session{iRegion,1},...
            xcorr_per_neuron_session{iRegion,2},4, 'filled')
        %makepretty_lite;
        %ylim([0.65, 0.8])
        makepretty;
        xlim([-0.2, 0.8])
        ylim([-0.2, 0.8])
        xlabel('correlation stim 1 v 2')
        ylabel('correlation stim 1 v 3')
        axis equal;
    %end
end



%% f.r. (each neuron)


%% neural trajectory 
% remove mean from each neuron (center the data) 
% (smooth) 
% PCA (3T * N matrix) 
% plot on PC1, 2, 3

for iRegion = 1:size(plot_regions,2)
    av_centered_smoothed_all_sessions= [];
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' ==1 | task_data_here.unitType' ==2) & ...
        task_data_here.pvalue_shuffled_005{1,2}' == 1;
        unitCount= 0;
   
    for iSession = 1:size(nonZero_idx,2)-1
        thisSession = nonZero_idx(iSession+1) -1;
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
        average_per_neuron = task_data_here.av_psth{2,thisSession}(these_units_session,:,:);
        
        if size(for_baseline_per_neuron_per_condition,2) == 39
            cond_inds = [10,34; 16,34; 16, 10];
        elseif size(for_baseline_per_neuron_per_condition,2) == 26
            cond_inds = [7,23; 11,23; 11, 7];
        else
            disp('wtf')
        end
        

        % center data (zscore)
        av1_centered = squeeze(average_per_neuron(:, cond_inds(1,1), :)) - nanmean(average_per_neuron(:, cond_inds(1,1), 1:200), 3)./...
            nanstd(average_per_neuron(:, cond_inds(1,1), 1:200), [], 3);
        av2_centered = squeeze(average_per_neuron(:, cond_inds(1,2), :)) - nanmean(average_per_neuron(:, cond_inds(1,2), 1:200), 3)./...
            nanstd(average_per_neuron(:, cond_inds(1,2), 1:200), [], 3);
        av3_centered = squeeze(average_per_neuron(:, cond_inds(2,1), :)) - nanmean(average_per_neuron(:, cond_inds(2,1), 1:200), 3)./...
            nanstd(average_per_neuron(:, cond_inds(2,1), 1:200), [], 3);

        % (smooth)
        av1_centered_smoothed = smoothdata(av1_centered, 2,'gaussian', [0, 0]);
        av2_centered_smoothed = smoothdata(av2_centered, 2,'gaussian', [0, 0]);
        av3_centered_smoothed = smoothdata(av3_centered, 2,'gaussian', [0, 0]);

        % matrix 
        if size(av1_centered_smoothed,2) == 1
            av_centered_smoothed_all = [av1_centered_smoothed', av2_centered_smoothed', av3_centered_smoothed'];
        else
            av_centered_smoothed_all = [av1_centered_smoothed, av2_centered_smoothed, av3_centered_smoothed];
        end
        % remove any nans 
        av_centered_smoothed_all(sum(isnan(av_centered_smoothed_all), 2) > 0, :) = [];

        av_centered_smoothed_all_sessions = [av_centered_smoothed_all_sessions; av_centered_smoothed_all];
    end
        % PCA 
        [pca_av,~,~,~,varex] = pca(av_centered_smoothed_all_sessions);



         figure(100); 
        subplot(3,3,(iRegion-1)*3+1)
        plot(-200:599, smoothdata(pca_av(1:800,1) - nanmean(pca_av(1:200,1)),'gaussian', [0 50]) ); hold on;
        plot(-200:599,  smoothdata(pca_av(801:1600,1) - nanmean(pca_av(801:1000,1)),'gaussian', [0 50]) );
        plot(-200:599,  smoothdata(pca_av(1601:2400,1) - nanmean(pca_av(1601:1800,1)),'gaussian', [0 50]) );
        ylabel(['PC1, expl. var.:', num2str(varex(1))])

        subplot(3,3,(iRegion-1)*3+2)
       plot(-200:599, smoothdata(pca_av(1:800,2) - nanmean(pca_av(1:200,2)),'gaussian', [0 50]) ); hold on;
        plot(-200:599,  smoothdata(pca_av(801:1600,2) - nanmean(pca_av(801:1000,2)),'gaussian', [0 50]) );
        plot(-200:599,  smoothdata(pca_av(1601:2400,2) - nanmean(pca_av(1601:1800,2)),'gaussian', [0 50]) );
        ylabel(['PC2, expl. var.:', num2str(varex(2))])

        subplot(3,3,(iRegion-1)*3+3)
         plot(-200:599, smoothdata(pca_av(1:800,3) - nanmean(pca_av(1:200,3)),'gaussian', [0 50]) ); hold on;
        plot(-200:599,  smoothdata(pca_av(801:1600,3) - nanmean(pca_av(801:1000,3)),'gaussian', [0 50]) );
        plot(-200:599,  smoothdata(pca_av(1601:2400,3) - nanmean(pca_av(1601:1800,3)),'gaussian', [0 50]) );
       ylabel(['PC3, expl. var.:', num2str(varex(3))])
       xlabel('time from stim (ms)')
       legend({'stim 1 (go)', 'stim 2 (go)', 'stim 3 (no-go)'})

       figure(102);
       subplot(1,3, iRegion); hold on;
       smBk =200;
       smFw=10;
       
       sm1 = smoothdata(pca_av(1:800,1) - nanmean(pca_av(1:200,1)),'gaussian', [smBk smFw]);
       sm2 = smoothdata(pca_av(1:800,2) - nanmean(pca_av(1:200,2)),'gaussian', [smBk smFw]) ;
       sm3 = smoothdata(pca_av(1:800,3) - nanmean(pca_av(1:200,3)),'gaussian', [smBk smFw]) ;

       sm1_2 = smoothdata(pca_av(801:1600,1) - nanmean(pca_av(801:1000,1)),'gaussian', [smBk smFw]);
       sm2_2 = smoothdata(pca_av(801:1600,2) - nanmean(pca_av(801:1000,2)),'gaussian', [smBk smFw]) ;
       sm3_2 = smoothdata(pca_av(801:1600,3) - nanmean(pca_av(801:1000,3)),'gaussian', [smBk smFw]) ;

       sm1_3 = smoothdata(pca_av(1601:2400,1) - nanmean(pca_av(1601:1800,1)),'gaussian', [smBk smFw]);
       sm2_3 = smoothdata(pca_av(1601:2400,2) - nanmean(pca_av(1601:1800,2)),'gaussian', [smBk smFw]) ;
       sm3_3 = smoothdata(pca_av(1601:2400,3) - nanmean(pca_av(1601:1800,3)),'gaussian', [smBk smFw]) ;

           plot3(sm1(200:end),sm2(200:end),sm3(200:end)); hold on;
           plot3(sm1_2(200:end),sm2_2(200:end),sm3_2(200:end));
           plot3(sm1_3(200:end),sm2_3(200:end),sm3_3(200:end))
           szMk=50
           cols = [rgb('Purple'); rgb('Green'); rgb('Orange');rgb('Black')];
           times =  [1,200, 700, 799];
           for idx=2:4
               iTime = times(idx);
               scatter3(sm1(iTime),sm2(iTime),sm3(iTime),szMk,cols(idx,:), 'filled')
               scatter3(sm1_2(iTime),sm2_2(iTime),sm3_2(iTime),szMk,cols(idx,:), 'filled')
               scatter3(sm1_3(iTime),sm2_3(iTime),sm3_3(iTime),szMk,cols(idx,:), 'filled')
           end
       

end


