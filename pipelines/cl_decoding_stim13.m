function [decodingAccuracy, correct_value, guessed_value] = cl_decoding_stim12(task_data, idx, keepVis, keepUnits, decoding_method, regularizationMethod, n_folds, plot_regions, plotMe)


%plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data.av_per_trial, 2)
    session_nCells(iSession) = size(task_data.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
trialTypes = [4, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);


task_data.animal_day_site_shank(isnan(task_data.animal_day_site_shank)) = 0;

thesePairs = [1, 2];
for iRegion = 1:size(plot_regions, 2)
    decoding_matrix = cell(size(nonZero_idx, 2)-1, 1);
    trialType_matrix = cell(size(nonZero_idx, 2)-1, 1);
    trialType = [];
    trialCount = 0;
    d_prime= [];

    if keepVis
        these_units = task_data.unit_area == iRegion & ...
            ismember(task_data.unitType', keepUnits) & ...
            task_data.pvalue_shuffled_005{1, idx}' == 1; %any(task_data.pvalue{1,idx} < 0.01,2);
    else
        these_units = task_data.unit_area == iRegion & ...
            ismember(task_data.unitType', keepUnits); % & ...
        %task_data.pvalue_shuffled_005{1,idx}' == 1;
    end
    unitCount = 0;
    for iSession = 1:size(nonZero_idx, 2) - 1
        thisSession = nonZero_idx(iSession+1) - 1;
        trials_no_move = task_data.trial_types{idx, thisSession}(task_data.no_move_trials{idx, thisSession}, :);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));


        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move, trialTypes, 'rows'); %QQ
            theseTrialTypes = trials_no_move(theseTrials, :);
            these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
            these_units_session_all = find(these_units);
            these_units_session_all = these_units_session_all( ...
                these_units_session_all >= session_cumCells(thisSession) & these_units_session_all <= session_cumCells(thisSession+1)-1);
            %figure();
            %i=2;
            %imagesc(squeeze(task_data.av_psth{idx,thisSession}(:,i,:)))
            activity_per_trial_neuron = task_data.av_per_trial{idx, thisSession}(these_units_session, theseTrials)';
            av_psth_here = task_data.av_psth{idx, thisSession}(these_units_session, :, :);
            for_baseline_per_neuron_per_condition = task_data.av_psth{idx, thisSession}(these_units_session, :, 1:200);
            if size(for_baseline_per_neuron_per_condition, 2) == 39
                cond_inds = [10, 16];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 26
                cond_inds = [7, 11];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 4
                cond_inds = [2, 3];
            else
                disp('wtf')
            end


            clearvars d_prime_session cond_fr keepIdx

            if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron > 0)))


                if ~isempty(activity_per_trial_neuron)
                    if size(activity_per_trial_neuron, 2) > 2 && size(activity_per_trial_neuron, 1) > 2

                        decoding_matrix{iSession} = activity_per_trial_neuron; % trials * neurons
                        trialType = [trialType; theseTrialTypes(:, 1)];
                        trialType_matrix{iSession} = theseTrialTypes(:, 1);
                         for iNeuron = 1:size(activity_per_trial_neuron, 2)
                        clearvars dprime_here
                        for iPair = 1
                            trials_1 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair), :), 'rows');
                            trials_2 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair), :), 'rows');

 average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ?
                                average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001));

                               
                        sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001); %./sqrt(sum(trials_1));
                        sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001); %./sqrt(sum(trials_2));
                        n_stim_1 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                        n_stim_2 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                        pooled_sd = sqrt((n_stim_1*sd_stim_1 * sd_stim_1 + n_stim_2*sd_stim_2 * sd_stim_2)./(n_stim_1+n_stim_2)) + 0.01;
                        dprime_here(:) = abs(average_stim_1-average_stim_2) / pooled_sd;
                        
                        end
                        d_prime = [d_prime, any(abs(dprime_here(:))>0.35)] ;
                         end
                               
                        %trialCount = trialCount + size(activity_per_trial_neuron,1);
                    end

                end
                unitCount = unitCount + size(activity_per_trial_neuron, 2);
            else
            end
        end

    end

    % combine all

    % Initialize total neuron and trial counts
    total_neurons = 0;
    total_trials = 1000;
    total_trials_1 = 1000;
    total_trials_2 = 1000;
    total_trials_3 = 1000;

    % Determine the total size
    for i = 1:numel(decoding_matrix)
        total_neurons = total_neurons + size(decoding_matrix{i}, 2);
       if ~isempty(decoding_matrix{i})
            total_trials = min(total_trials, size(decoding_matrix{i}, 1));
            total_trials_1 =  min(total_trials_1, sum(trialType_matrix{i} == 4));
            total_trials_2 =  min(total_trials_2, sum(trialType_matrix{i} == 6));
           
       end
    end

    % Initialize the merged matrix with NaNs
    mergedMatrix = NaN(total_trials, total_neurons);

    % Variables to keep track of the current position in the matrix
    current_neuron_index = 1;
    current_trial_index = 1;

    % Fill the merged matrix
    for i = 1:numel(decoding_matrix)
        currentMatrix = decoding_matrix{i};
        n_neurons = size(currentMatrix, 2);
        n_trials = size(currentMatrix, 1);

        % Place the current matrix in the merged matrix
        if ~isempty(currentMatrix)
            currentMatrix_thisStim = currentMatrix(trialType_matrix{i} == 4, :);
            mergedMatrix(1:total_trials_1, current_neuron_index:(current_neuron_index + n_neurons - 1)) = currentMatrix_thisStim(1:total_trials_1, :);
            
            currentMatrix_thisStim = currentMatrix(trialType_matrix{i} == 6, :);
            mergedMatrix(total_trials_1+1:total_trials_1+total_trials_2, current_neuron_index:(current_neuron_index + n_neurons - 1)) = currentMatrix_thisStim(1:total_trials_2, :);
            
          end

        % Update the indices
        current_neuron_index = current_neuron_index + n_neurons;
        current_trial_index = current_trial_index + n_trials;
    end
    trialTypes_merged = [ones(total_trials_1,1); ones(total_trials_2,1).*2];

    % % TEST: fit the whole thing 
    % fit3 = glmnet(mergedMatrix, trialTypes_merged, 'multinomial');
    % glmnetPrint(fit3);
    % prediction_cv = glmnetPredict(fit3, mergedMatrix, [], 'class'); %extract coefficients at a single value of lambda%    
    % prediction_cv_highlambda = prediction_cv(:,100);
    % correct = trialTypes_merged(prediction_cv_highlambda == trialTypes_merged);
    % correct_perStim(1) = sum(correct==1)/total_trials_1;
    % correct_perStim(2) = sum(correct==2)/total_trials_2;
    % correct_perStim(3) = sum(correct==3)/total_trials_3;

    % real thing: divide into 5 folds. 
    %n_folds = 20;
    fold_length_1 = floor(total_trials_1/n_folds);
    fold_length_2 = floor(total_trials_2/n_folds);
   
    folds_trials = [];
    for i = 1:2
        fold_length = eval(['fold_length_' num2str(i)]);
        total_trials = eval(['total_trials_' num2str(i)]);
        for j = 1:n_folds-1
            folds_trials = [folds_trials; ones(fold_length, 1) * j];
        end
        folds_trials = [folds_trials; ones(total_trials - fold_length * (n_folds-1), 1) * n_folds];
    end

    % compute for each fold on 4/5th, 1/5 
    correct_perStim_cv = nan(2, n_folds);
    for iFold = 1:n_folds
        trainSet = mergedMatrix(folds_trials ~= iFold, :);
        train_trials = trialTypes_merged(folds_trials ~= iFold, :);
        testSet = mergedMatrix(folds_trials == iFold, :);
        test_trials = trialTypes_merged(folds_trials == iFold, :);
        
        if strcmp(regularizationMethod, 'L1')
            cvfit = cvglmnet(trainSet, train_trials, 'multinomial');
    
            % Best lambda
            lambda_optimal = cvfit.lambda_min;
    
            options = glmnetSet;
            options.lambda = lambda_optimal;
            fit3 = glmnet(trainSet, train_trials, 'multinomial', options);
            %glmnetPrint(fit3);
            prediction_cv = glmnetPredict(fit3, testSet, lambda_optimal, 'class'); %extract coefficients at a single value of lambda%

            ncoef = glmnetCoef(fit3,[],false);
           % figure(); hold on;
           % plot(ncoef{1});plot(ncoef{2});plot(ncoef{3});
        else
           % keepNeurons = d_prime;
            options = glmnetSet;
            fit3 = glmnet(trainSet(:,logical(d_prime)), train_trials, 'multinomial', options);
            lambda_optimal = 0;
            prediction_cv = glmnetPredict(fit3, testSet(:, logical(d_prime)), lambda_optimal, 'class'); %extract coefficients at a single value of lambda% 
    
            
        end

        %prediction_cv_highlambda = prediction_cv(:,100);

        correct = test_trials(prediction_cv == test_trials);
        correct_perStim_cv(1,iFold) = sum(correct==1)/sum(test_trials==1);
        correct_perStim_cv(2,iFold) = sum(correct==2)/sum(test_trials==2);
        correct_value{iRegion}{iFold} = test_trials;
        guessed_value{iRegion}{iFold} = prediction_cv;

    end
    
  decodingAccuracy(iRegion, :, :) =  correct_perStim_cv;
  %correct_value{iFold} = test_trials;
  %guessed_value{iFold} = prediction_cv;

end

