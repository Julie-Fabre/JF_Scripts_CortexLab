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

%task_data_here = load('/home/julie/task_data_goNogo2.mat');

plot_regions = [1, 2, 5]; % Striatum, GPe, SNr 
for iSession = 1:size(task_data_here.av_per_trial,2)
    session_nCells(iSession) = size(task_data_here.av_per_trial{2,iSession},1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells>0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
% session_cumCells = cumsum(session_nCells);

% glmnet (logistic regression, L1 penalty, cross validated)
for iRegion = 1:size(plot_regions,3)
    %these_units = task_data_here.unit_area == iRegion & ...
    %    (task_data_here.unitType' ==1 | task_data_here.unitType' ==2);
   
    for iSession = 1:size(nonZero_idx,2)
        thisSession = nonZero_idx(iSession);
        trials_no_move = task_data_here.trial_types{2,thisSession}(task_data_here.no_move_trials{2,thisSession},:);
        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move,...
                trialTypes, 'rows');%QQ
            theseTrialTypes = trials_no_move(theseTrials,:);
           % if session_nCells > 0 & 
            %these_units_session = these_units(session_nCells(thisSession):session_nCells(thisSession+1));
            %if ~isempty(these_units_session)
                activity_per_trial_neuron = task_data_here.av_per_trial{2,thisSession}(:,theseTrials)';
        %cvglmnet(activity_per_trial_neuron, trial_type, family, options  type,
        %nfolds)
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
           % end
        end

    end
    %g2=randsample(2,100,true);
    %x=randn(100,20);
    %fit2=glmnet(x,g2,'binomial');
end

% glmnet (logistic regression, selecting top x% cells, cross validated)
fit3=glmnet(activity_per_trial_neuron(), theseTrialTypes(:,1), 'multinomial');

x=randn(100,20);
g4=randsample(4,100,true);
fit3=glmnet(x,g4,'multinomial');

opts=struct('mtype','grouped');
fit3a=glmnet(x,g4,'multinomial',opts);

% dot average product activity of each cell across trials for each time point


   n=500; p=30;
   nzc=fix(p/10);
   x=randn(n,p);
   beta3=randn(10,3);
   beta3=cat(1,beta3,zeros(p-10,3));
   f3=x*beta3;
   p3=exp(f3);
   p3=bsxfun(@rdivide,p3,sum(p3,2));
   g3=mnrnd(1,p3);
   g3=g3*(1:size(p3,2))';
   cvfit=cvglmnet(x,g3,'multinomial', [], [] ,5);
   cvglmnetPlot(cvfit);


      x=randn(100,20);
   y=randn(100,1);
   fit1 = glmnet(x,y);
   glmnetPrint(fit1);
   glmnetPredict(fit1,[],0.01,'coef')  %extract coefficients at a single value of lambda
   glmnetPredict(fit1,x(1:10,:),[0.01,0.005]') 


    g4=randsample(4,100,true);
    fit3=glmnet(x,g4,'multinomial');
    opts=struct('mtype','grouped');
    fit3a=glmnet(x,g4,'multinomial',opts);
       glmnetPrint(fit3);
   glmnetPredict(fit3,[],0.01,'coef')  %extract coefficients at a single value of lambda
   glmnetPredict(fit3,x(1:10,:),[0.01,0.005]') 

