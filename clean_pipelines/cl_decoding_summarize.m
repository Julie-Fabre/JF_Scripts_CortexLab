%% cl_decoding_summarize

% parameters
keepVis = 1;
keepUnits = [1, 2]; % 1:good, 2:mua, 3:non-somatic, 4:noise
plotMe = false;
plot_regions = [1, 2, 3];

% initialize variables
d_prime = cell(3, 1);
d_prime_session_num = cell(3, 1);
%d_prime_animal_num = cell(3, 1);
d_prime_session_fraction = cell(3, 1);
d_prime_session_faction_median = cell(3, 1);

% run loop
n_folds = 10;
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive.mat');
        idx = 5;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
        idx = 2;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
        idx = 2;

    end

    [decodingAccuracy{iTask}, correctValue{iTask}, guessedValue{iTask}] ...
        = cl_decoding(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'L1', n_folds, plot_regions, plotMe);

    [decodingAccuracy12{iTask}, correctValue12{iTask}, guessedValue12{iTask}] ...
        = cl_decoding_stim12(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'L1', n_folds, plot_regions, plotMe);

      [decodingAccuracy13{iTask}, correctValue13{iTask}, guessedValue13{iTask}] ...
        = cl_decoding_stim13(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'L1', n_folds, plot_regions, plotMe);

        [decodingAccuracy23{iTask}, correctValue23{iTask}, guessedValue23{iTask}] ...
        = cl_decoding_stim23(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'L1', n_folds, plot_regions, plotMe);

end

%% all stims
% comparing tasks side by side
figure();
for iRegion =1:3
for iStim=1:3 
    subplot(3, 3, iStim+(iRegion - 1)*(3));
    violinplot([squeeze(decodingAccuracy{1}(iRegion, iStim, :)), squeeze(decodingAccuracy{2}(iRegion, iStim, :)),...
        squeeze(decodingAccuracy{3}(iRegion, iStim, :))], [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 
    title([num2str(iStim), ', ', num2str(iRegion)])
    xlabel('Task #')
    ylabel(['region #' num2str(iRegion)])
end
end

% comparing regions side by side 
figure();
for iTask = 1:3
for iStim=1:3 
    subplot(3, 3, iStim+(iTask- 1)*(3));
    violinplot([squeeze(decodingAccuracy{iTask}(1, iStim, :)), squeeze(decodingAccuracy{iTask}(2, iStim, :)),...
        squeeze(decodingAccuracy{iTask}(3, iStim, :))], [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 
    title([num2str(iStim), ', ', num2str(iTask)])
    xlabel('Region #')
    ylabel(['task #' num2str(iTask)])
end
end

% decoding matrix 
figure();
for iTask =1:3
for iRegion =1:3
    clearvars  matrixDecoding
    %subplot(3, 3, iStim+(iRegion - 1)*(3));
    for iFold = 1:10
        matrixDecoding(iFold, 1, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 1, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 1, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);

        matrixDecoding(iFold, 2, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 2, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 2, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);

        matrixDecoding(iFold, 3, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 3, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 3, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
    end
    matrixDecoding_simple = squeeze(nanmean(matrixDecoding, 1));
    subplot(3, 3, iTask+(iRegion - 1)*(3));
    imagesc(matrixDecoding_simple)
    colormap(brewermap([], 'YlOrRd'));
    clim([0 1])
end
end
prettify_plot;

%% stims 1 and 2 



for iTask =1:3

figure();
%for iStim=1:2
   % subplot(3, 3, iStim+(iTask - 1)*(3));
    violinplot([squeeze(nanmean(decodingAccuracy12{iTask}(1, :, :))), ...
        squeeze(nanmean(decodingAccuracy12{iTask}(2, :, :))), ....
        squeeze(nanmean(decodingAccuracy12{iTask}(3, :, :)))]...
        , [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 

    
     [p, t, stats] = anova1([squeeze(nanmean(decodingAccuracy12{iTask}(1, :, :))); ...
         squeeze(nanmean(decodingAccuracy12{iTask}(2, :, :))); ....
         squeeze(nanmean(decodingAccuracy12{iTask}(3, :, :)))], [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]);
[c, m, h, gnames] = multcompare(stats);
   
    % stim12_table = table;
    % stim12_table.decoding = [squeeze(nanmean(decodingAccuracy12{iTask}(1, :, :))); ...
    %     squeeze(nanmean(decodingAccuracy12{iTask}(2, :, :))); ....
    %     squeeze(nanmean(decodingAccuracy12{iTask}(3, :, :)))];
    % stim12_table.region = [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]; 
    % lme(stim12_table, 'decoding ~ 1 + region')
    
%end
end

for iTask =1:3

figure();
%for iStim=1:2
   % subplot(3, 3, iStim+(iTask - 1)*(3));
    violinplot([squeeze(nanmean(decodingAccuracy13{iTask}(1, :, :))), ...
        squeeze(nanmean(decodingAccuracy13{iTask}(2, :, :))), ....
        squeeze(nanmean(decodingAccuracy13{iTask}(3, :, :)))]...
        , [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 

    
     [p, t, stats] = anova1([squeeze(nanmean(decodingAccuracy13{iTask}(1, :, :))); ...
         squeeze(nanmean(decodingAccuracy13{iTask}(2, :, :))); ....
         squeeze(nanmean(decodingAccuracy13{iTask}(3, :, :)))], [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]);
[c, m, h, gnames] = multcompare(stats);
   
    % stim12_table = table;
    % stim12_table.decoding = [squeeze(nanmean(decodingAccuracy12{iTask}(1, :, :))); ...
    %     squeeze(nanmean(decodingAccuracy12{iTask}(2, :, :))); ....
    %     squeeze(nanmean(decodingAccuracy12{iTask}(3, :, :)))];
    % stim12_table.region = [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]; 
    % lme(stim12_table, 'decoding ~ 1 + region')
    
%end
end


for iTask =1:3

figure();
%for iStim=1:2
   % subplot(3, 3, iStim+(iTask - 1)*(3));
    violinplot([squeeze(nanmean(decodingAccuracy23{iTask}(1, :, :))), ...
        squeeze(nanmean(decodingAccuracy23{iTask}(2, :, :))), ....
        squeeze(nanmean(decodingAccuracy23{iTask}(3, :, :)))]...
        , [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 

    
     [p, t, stats] = anova1([squeeze(nanmean(decodingAccuracy23{iTask}(1, :, :))); ...
         squeeze(nanmean(decodingAccuracy23{iTask}(2, :, :))); ....
         squeeze(nanmean(decodingAccuracy23{iTask}(3, :, :)))], [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]);
[c, m, h, gnames] = multcompare(stats);
   
    % stim12_table = table;
    % stim12_table.decoding = [squeeze(nanmean(decodingAccuracy12{iTask}(1, :, :))); ...
    %     squeeze(nanmean(decodingAccuracy12{iTask}(2, :, :))); ....
    %     squeeze(nanmean(decodingAccuracy12{iTask}(3, :, :)))];
    % stim12_table.region = [ones(n_folds,1); ones(n_folds,1).*2; ones(n_folds,1).*3]; 
    % lme(stim12_table, 'decoding ~ 1 + region')
    
%end
end



% decoding matrix 
figure();
for iTask =1:3
for iRegion =1:3
    clearvars  matrixDecoding
    %subplot(3, 3, iStim+(iRegion - 1)*(3));
    for iFold = 1:10
        matrixDecoding(iFold, 1, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 1, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 1, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 1 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);

        matrixDecoding(iFold, 2, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 2, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 2, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 2 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);

        matrixDecoding(iFold, 3, 1) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 1)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 3, 2) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 2)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
        matrixDecoding(iFold, 3, 3) = sum(correctValue{iTask}{iRegion}{iFold} == 3 & guessedValue{iTask}{iRegion}{iFold} == 3)./...
            sum(correctValue{iTask}{iRegion}{iFold} == 1);
    end
    matrixDecoding_simple = squeeze(nanmean(matrixDecoding, 1));
    subplot(3, 3, iTask+(iRegion - 1)*(3));
    imagesc(matrixDecoding_simple)
    colormap(brewermap([], 'YlOrRd'));
    clim([0 1])
end
end
prettify_plot;