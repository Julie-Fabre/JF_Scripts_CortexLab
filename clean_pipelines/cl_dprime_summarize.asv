
%% cl_dprime_summarize

% parameters
keepVis = 1;
keepUnits = [1, 2]; % 1:good, 2:mua, 3:non-somatic, 4:noise
plotMe = false;
plotRegions = [1,2,3];

% initialize variables
d_prime = cell(3, 1);
d_prime_session_num = cell(3, 1);
d_prime_animal_num = cell(3, 1);
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

    [d_prime{iTask}, d_prime_session_num{iTask}, d_prime_animal_num{iTask}, d_prime_session_fraction{iTask}]...
        = cl_dprime(task_data, idx, keepVis, keepUnits, plotRegions, plotMe);
end

% plot overlaid histograms
cols = ya_getColors(3);
figure();
for iRegion = 1:3
    for iPair = 1:3
        for iTask = 3:-1:1
            subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)));
            hold on;
            kp = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, iPair))));

            histogram(d_prime{iTask}{iRegion}(kp, iPair), 0:0.05:2, 'Normalization', 'probability',...
                'FaceColor', cols(iTask,:), 'EdgeColor', cols(iTask,:), 'FaceAlpha', 0.1);
            xlabel('d-prime')
            ylabel('# of neurons')
            xlim([0, 2])
        end
    end
end
legend({'gonogo', 'gogogo', 'naive'})
prettify_plot('YLimits', 'all')

% plot histograms sorted by session and mouse (apris next to each other)

% plot histograms


% plot