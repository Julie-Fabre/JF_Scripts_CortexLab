%% cl_decoding_summarize

% parameters
keepVis = 1;
keepUnits = [1, 2]; % 1:good, 2:mua, 3:non-somatic, 4:noise
plotMe = false;
plotRegions = [1, 2, 3];

% initialize variables
d_prime = cell(3, 1);
d_prime_session_num = cell(3, 1);
%d_prime_animal_num = cell(3, 1);
d_prime_session_fraction = cell(3, 1);
d_prime_session_faction_median = cell(3, 1);

% run loop
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

    [decodingAccuracy{iTask}] ...
        = cl_decoding(task_data, idx, keepVis, keepUnits, 'cvglmnet', 'elastic', plot_regions, plotMe);

end

figure();
for iRegion =1:3
for iStim=1:3 
    subplot(3, 3, iStim+(iRegion - 1)*(3));
    scatter(ones(5,1), squeeze(decodingAccuracy{1}(iRegion, iStim, :))); 
    hold on;
    scatter(ones(5,1).*2, squeeze(decodingAccuracy{2}(iRegion, iStim, :))); 
    hold on;
    scatter(ones(5,1).*3, squeeze(decodingAccuracy{3}(iRegion, iStim, :))); 
end
end