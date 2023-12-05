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
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo3.mat');
        idx = 2;

    end

    [decodingAccuracy{iTask}, correctValue{iTask}, guessedValue{iTask}] ...
        = cl_decoding(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'L1', n_folds, plot_regions, plotMe);

end

figure();
for iRegion =1:3
for iStim=1:3 
    subplot(3, 3, iStim+(iRegion - 1)*(3));
    violinplot([squeeze(decodingAccuracy{1}(iRegion, iStim, :)), squeeze(decodingAccuracy{2}(iRegion, iStim, :)),...
        squeeze(decodingAccuracy{3}(iRegion, iStim, :))], [ones(n_folds,1), ones(n_folds,1).*2, ones(n_folds,1).*3]); 
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