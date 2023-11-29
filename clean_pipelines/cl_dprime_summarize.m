
%% cl_dprime_summarize

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

    [d_prime{iTask}, d_prime_session_num{iTask}, d_prime_animal_num{iTask}, d_prime_session_fraction{iTask}] ...
        = cl_dprime(task_data, idx, keepVis, keepUnits, plotRegions, plotMe);
end

% plot overlaid histograms
cols = ya_getColors(3);
figure();
for iRegion = 1:3
    for iPair = 1:3
        for iTask = 3:-1:1
            if iTask == 1 && (iRegion == 2 || iRegion == 3)
                continue;
            end
            subplot(3, size(plotRegions, 2), iPair+(iRegion - 1)*(size(plotRegions, 2)));
            hold on;
            kp = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, iPair))));

            histogram(d_prime{iTask}{iRegion}(kp, iPair), 0:0.05:4, 'Normalization', 'probability', ...
                'FaceColor', cols(iTask, :), 'EdgeColor', cols(iTask, :), 'FaceAlpha', 0.1);
            xlabel('d-prime')
            ylabel('# of neurons')
            xlim([0, 2])
            if iTask == 1 && iRegion == 1 && iPair == 1
                legend({'gonogo', 'gogogo', 'naive'})
            end
        end
    end
end
prettify_plot('YLimits', 'all')

% plot histograms sorted by session and mouse (pairs next to each other)
figure();
for iPair = 1:3
    for iRegion = 1:3
        subplot(3, 3, iPair+(iRegion - 1)*3);
        num_recs = size(unique(d_prime_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];


        iTask = 2;
        kp1 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));
        nums1 = d_prime_session_num{iTask}{iRegion}(kp1, 1);
        [~, ~, ic] = unique(nums1, 'stable');
        nums_u1 = ic';

        iTask = 3;
        kp2 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));
        nums2 = d_prime_session_num{iTask}{iRegion}(kp2, 1);
        [~, ~, ic] = unique(nums2, 'stable');
        nums_u2 = ic';

        theseColors = repelem(cols(2, :), length(unique(nums_u1)), 1);
        theseColors = [theseColors; repelem(cols(3, :), length(unique(nums_u2)), 1)];

        violinplot([d_prime{2}{iRegion}(kp1, iPair); d_prime{3}{iRegion}(kp2, iPair)], ...
            [nums_u1'; max(nums_u1) + nums_u2'], ...
            'ViolinColor', theseColors);
        ylim([0, 2.6])

        uu = unique(d_prime_session_num{2}{iRegion});
        mouseys = [d_prime_animal_num{2}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{2}{iRegion}(:, 1))), 1) ...
            }];
        mouseys_u = unique(mouseys);
        colsMice = ya_getColors(max(mouseys_u));
        cc = 0;
        for iMousey = 1:size(mouseys_u, 2)
            idx_m = find(mouseys == mouseys_u(iMousey));
            line([cc, cc + length(idx_m) + 0.5], [4, 4], 'Color', colsMice(mouseys_u(iMousey), :));
            cc = cc + length(idx_m);
        end

        uu = unique(d_prime_session_num{3}{iRegion});
        mouseys = [d_prime_animal_num{3}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{3}{iRegion}(:, 1))), 1) ...
            }];

        mouseys_u = unique(mouseys);
        colsMice = ya_getColors(max(mouseys_u));
        %cc =0;
        for iMousey = 1:size(mouseys_u, 2)
            idx_m = find(mouseys == mouseys_u(iMousey));
            line([cc, cc + length(idx_m) + 0.5], [4, 4], 'Color', colsMice(mouseys_u(iMousey), :));
            cc = cc + length(idx_m);
        end

        axis tight;


    end
end
prettify_plot('YLimits', 'all')


% general
theseColors = ya_getColors(3);
figure();
for iPair = 1:3
    for iRegion = 1:3
        subplot(3, 3, iPair+(iRegion - 1)*3);
        num_recs = size(unique(d_prime_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];

        iTask = 1;
        kp0 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 2;
        kp1 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 3;
        kp2 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        violinplot([d_prime{1}{iRegion}(kp0, iPair); ...
            d_prime{2}{iRegion}(kp1, iPair); d_prime{3}{iRegion}(kp2, iPair)], ...
            [ones(length(d_prime{1}{iRegion}(kp0, iPair)), 1); ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors);
        ylim([0, 4])

    end
end

%% ~~ THE TESTS ~~
%% NORMAL
theseColors = ya_getColors(3);
figure();
for iPair = 1:3
    for iRegion = 1:3
        subplot(3, 3, iPair+(iRegion - 1)*3);
        num_recs = size(unique(d_prime_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];

        iTask = 1;
        kp0 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 2;
        kp1 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 3;
        kp2 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        %
        % violinplot([sqrt(d_prime{1}{iRegion}(kp0, iPair));...
        %     sqrt(d_prime{2}{iRegion}(kp1, iPair)); d_prime{3}{iRegion}(kp2, iPair)], ...
        %     [ones(length(d_prime{1}{iRegion}(kp0, iPair)),1);...
        %     ones(length(d_prime{2}{iRegion}(kp1, iPair)),1).*2; ...
        %     ones(length(d_prime{3}{iRegion}(kp2, iPair)),1).*3], ...
        %     'ViolinColor', theseColors);
        violinplot([ ...
            (d_prime{2}{iRegion}(kp1, iPair)); (d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors);

        lmeTable = table;
        lmeTable.dprime = (d_prime{2}{iRegion}(kp1, iPair));
        lmeTable.sessionType = ones(size(d_prime{2}{iRegion}(kp1, iPair), 1), 1);
        lmeTable.session = d_prime_session_num{2}{iRegion}(kp1, iPair);

        uu = unique(d_prime_session_num{2}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{2}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{2}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{2}{iRegion}(:, iPair))), 1) ...
            }];

        mousee = d_prime_session_num{2}{iRegion}(kp1, iPair);
        uu = unique(d_prime_session_num{2}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = mousee;
        
        lmeTable2 = table;

        uu = unique(d_prime_session_num{3}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{3}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{3}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{3}{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable2.dprime = (d_prime{3}{iRegion}(kp2, iPair));
        lmeTable2.sessionType = ones(size(d_prime{3}{iRegion}(kp2, iPair), 1), 1).*2;
        lmeTable2.session = d_prime_session_num{3}{iRegion}(kp2, iPair);
        mousee = d_prime_session_num{3}{iRegion}(kp2, iPair);
        uu = unique(d_prime_session_num{3}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable2.mouse = mousee;

        lmeTable = [lmeTable;lmeTable2];
        lme = fitlme(lmeTable, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        %lme = fitlme(lmeTable, 'dprime ~ sessionType + (1|session)'); %what if I add 1|neuron?
       ylim([0, 4.2])
        prettify_pvalues(gca, [1], [2], lme.coefTest)

        

    end
end




%% SQRT 
theseColors = ya_getColors(3);
figure();
for iPair = 1:3
    for iRegion = 1:3
        subplot(3, 3, iPair+(iRegion - 1)*3);
        num_recs = size(unique(d_prime_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];

        iTask = 1;
        kp0 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 2;
        kp1 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ...
            ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        iTask = 3;
        kp2 = find(abs(d_prime{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{iRegion}(:, 1))) & ~isnan(abs(d_prime{iTask}{iRegion}(:, 1))));


        %
        % violinplot([sqrt(d_prime{1}{iRegion}(kp0, iPair));...
        %     sqrt(d_prime{2}{iRegion}(kp1, iPair)); d_prime{3}{iRegion}(kp2, iPair)], ...
        %     [ones(length(d_prime{1}{iRegion}(kp0, iPair)),1);...
        %     ones(length(d_prime{2}{iRegion}(kp1, iPair)),1).*2; ...
        %     ones(length(d_prime{3}{iRegion}(kp2, iPair)),1).*3], ...
        %     'ViolinColor', theseColors);
        violinplot([ ...
            sqrt(d_prime{2}{iRegion}(kp1, iPair)); sqrt(d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors);

        lmeTable = table;
        lmeTable.dprime = sqrt(d_prime{2}{iRegion}(kp1, iPair));
        lmeTable.sessionType = ones(size(d_prime{2}{iRegion}(kp1, iPair), 1), 1);
        lmeTable.session = d_prime_session_num{2}{iRegion}(kp1, iPair);

        uu = unique(d_prime_session_num{2}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{2}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{2}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{2}{iRegion}(:, iPair))), 1) ...
            }];

        mousee = d_prime_session_num{2}{iRegion}(kp1, iPair);
        uu = unique(d_prime_session_num{2}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = mousee;
        
        lmeTable2 = table;

        uu = unique(d_prime_session_num{3}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{3}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{3}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{3}{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable2.dprime = sqrt(d_prime{3}{iRegion}(kp2, iPair));
        lmeTable2.sessionType = ones(size(d_prime{3}{iRegion}(kp2, iPair), 1), 1).*2;
        lmeTable2.session = d_prime_session_num{3}{iRegion}(kp2, iPair);
        mousee = d_prime_session_num{3}{iRegion}(kp2, iPair);
        uu = unique(d_prime_session_num{3}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable2.mouse = mousee;

        lmeTable = [lmeTable;lmeTable2];
        lme = fitlme(lmeTable, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        %lme = fitlme(lmeTable, 'dprime ~ sessionType + (1|session)'); %what if I add 1|neuron?
       
        prettify_pvalues(gca, [1], [2], lme.coefTest)

        ylim([0, 2.1])

    end
end




% plot
% plot histograms sorted by session and mouse (tasks next to each other)

% plot histograms sorted by session and mouse (regions next to each other)
