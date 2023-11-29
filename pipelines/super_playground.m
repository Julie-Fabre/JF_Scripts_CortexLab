
%% super playground

JF_load_AP;

% for each unit below depth, get psth.
keepUnits = find(template_depths > 1500);

trial_conditions_clean = trial_conditions;
trial_conditions_clean(trial_conditions(:, 1) == 4, 1) = 100; %go1
trial_conditions_clean(trial_conditions(:, 1) == 12, 1) = 101; %go2
trial_conditions_clean(trial_conditions(:, 1) == 6, 1) = 102; %no go
trial_conditions_clean(trial_conditions(:, 1) == 11, 1) = 103; %no like
trial_conditions_clean(trial_conditions(:, 1) == 13, 1) = 104; %go like
theseImages = 100:102;
%trial_conditions(trial_conditions(:,1)>13,1) = trial_conditions(trial_conditions(:,1)>13,1)-13;

theseImages_trials = ismember(trial_conditions_clean(:, 1), theseImages) & ismember(trial_conditions_clean(:, 2), -90) & no_move_trials;

theseTimes = stimOn_times(theseImages_trials);

raster_window = [-0.5, 1];
psth_bin_size = 0.01;
for iKeepUnit = 1:length(keepUnits)

    iKeepUnit = iKeepUnit + 1;
    
    [curr_psth, curr_raster, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
        keepUnits(iKeepUnit), raster_window, psth_bin_size, ...
        theseTimes, trial_conditions_clean(theseImages_trials, 1));

    clf;
    plot(squeeze(curr_psth(1,:)))

    % average firing rate
    startIdx = find(t >= 0.05, 1, 'first');
    stopIdx = find(t >= 0.15, 1, 'first');
    av_per_trial(iKeepUnit, :) = nanmean(curr_raster(:, startIdx:stopIdx), 2);
    % psth
    psth_per_cell(iKeepUnit, :, :) = curr_psth;

end

figure();
subplot(131)
imagesc(squeeze(psth_per_cell(:, 1, :))); hold on;
title('stim 1 (go)')
xlabel('time from stim(s)')
ylabel('neuron #')
colorbar;
subplot(132)
imagesc(squeeze(psth_per_cell(:, 2, :))); hold on;
title('stim 2 (go)')
colorbar;
subplot(133)
imagesc(squeeze(psth_per_cell(:, 3, :))); hold on;
title('stim 3 (no go)')
colorbar;

figure();
imagesc(nanmean(av_per_trial,2))

family = 'multinomial'; %'binomial'; %logistic
% fit = glmnet(activity_per_trial_neuron, theseTrialTypes, family, options);
% get a good lambda value (?)

theseTrialTypes =  trial_conditions_clean(theseImages_trials, 1);
folds_trials = randsample(5, size(theseTrialTypes, 1), true);

 activity_per_trial_neuron = av_per_trial';
fit_cv_lambda = cvglmnet(activity_per_trial_neuron, ...
    theseTrialTypes(:, 1), family, [], [], 5);

use_lambda = fit_cv_lambda.lambda_min;
%fit_cv_lambda.glmnet_fit

%cvglmnetPlot(fit_cv_lambda);
%cvglmnetPredict(fit_cv_lambda,activity_per_trial_neuron(:,:),'lambda_min')

trialTypes = theseImages;
for iFold = 1:5 %divide into 5, train on 4/5th and test on remaining 1/5th

    theseFolds = ones(5, 1);
    theseFolds(iFold) = 0;
    options.lambda = use_lambda;
    fit_cv = glmnet(activity_per_trial_neuron(ismember(folds_trials, find(theseFolds)), :), ...
        theseTrialTypes(ismember(folds_trials, find(theseFolds)), 1), family, options);
    %glmnetPredict(object, newx, s, type, exact, offset)
    prediction_cv = glmnetPredict(fit_cv, activity_per_trial_neuron(ismember(folds_trials, iFold), :), [], [], true); %, use_lambda, 'link');
    [~ ,trials_pred_thisType] = max(prediction_cv');
    
    for iTrialType = 1:3
        trials_thisType = double(theseTrialTypes(ismember(folds_trials, iFold), 1) == trialTypes(iTrialType)) .* iTrialType;

        accuracy_cv(iFold, iTrialType) = sum(trials_thisType == trials_pred_thisType') ./ length(find(trials_thisType)); %which trials correctly labelled
    end


end

for ii = 100
figure();
plot(fit_cv_lambda.glmnet_fit.beta{1}(:,ii))
end
disp(accuracy_cv)
figure();

imagesc(accuracy_cv)
colorbar;
xlabel('stim type')
ylabel('validation fold')

figure();
scatter([1,1,1,1,1], accuracy_cv(:, 1), 'filled'); hold on;
boxplot( [accuracy_cv(:, 1), accuracy_cv(:, 2), accuracy_cv(:, 3)], [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
scatter([2,2,2,2,2], accuracy_cv(:, 2), 'filled'); hold on;
scatter([3,3,3,3,3], accuracy_cv(:, 3), 'filled'); hold on;

xlabel('stim type')
ylabel('fraction of trials correctly classified')
makepretty;
xlim([0.5, 3.5])

