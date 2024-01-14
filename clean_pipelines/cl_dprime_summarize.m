
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

% initialize variables
fr = cell(3, 1);
fr_session_num = cell(3, 1);
%d_prime_animal_num = cell(3, 1);
fr_session_fraction = cell(3, 1);
fr_session_faction_median = cell(3, 1);

% run loop
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat');
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

     [fr{iTask}, fr_session_num{iTask}, fr_animal_num{iTask}, fr_session_fraction{iTask}] ...
        = cl_frPerStim(task_data, idx, keepVis, keepUnits, plotRegions, plotMe);
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

        violinplot([d_prime{2}{iRegion}(kp1, iPair); d_prime{iTask}{iRegion}(kp2, iPair)], ...
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

        uu = unique(d_prime_session_num{iTask}{iRegion});
        mouseys = [d_prime_animal_num{iTask}{iRegion, ...
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
        if iRegion ==1
        violinplot([ ...
            sqrt(d_prime{1}{iRegion}(kp0, iPair));sqrt(d_prime{2}{iRegion}(kp1, iPair)); sqrt(d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ones(length(d_prime{1}{iRegion}(kp0, iPair)), 1) .* 1;
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([3,1,2],:));
        else
             violinplot([ ...
          sqrt(d_prime{2}{iRegion}(kp1, iPair)); sqrt(d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([1,2],:));
     
        end

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

        if iRegion ==1
        lmeTable3 = table;

        uu = unique(d_prime_session_num{1}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{1}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable3.dprime = sqrt(d_prime{1}{iRegion}(kp0, iPair));
        lmeTable3.sessionType = ones(size(d_prime{1}{iRegion}(kp0, iPair), 1), 1).*3;
        lmeTable3.session = d_prime_session_num{1}{iRegion}(kp0, iPair);
        mousee = d_prime_session_num{1}{iRegion}(kp0, iPair);
        uu = unique(d_prime_session_num{1}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable3.mouse = mousee;

        lmeTable13 = [lmeTable3;lmeTable2];
        lmeTable12 = [lmeTable3;lmeTable];
        end

        lmeTable23 = [lmeTable;lmeTable2];
        

        lme23 = fitlme(lmeTable23, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        if iRegion==1
        try
        lme13 = fitlme(lmeTable13, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        lme12 = fitlme(lmeTable12, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        
        catch
            lme13.coefTest=NaN;
            lme12.coefTest=NaN;
        end
        end
       
        
        %contrastMatrix = [0 1 0; 0 0 1]'; % For Type1 vs Type2, and Type1 vs Type3 comparisons
        %pValue = coefTest(lme, contrastMatrix);
        
        ylim([-0.5, 2.5])
        if iRegion==1
        prettify_pvalues(gca, [1,2,1], [2,3,3], [lme12.coefTest, lme23.coefTest, lme13.coefTest])
        else
             prettify_pvalues(gca, [1], [2], [lme23.coefTest])
       
        end


    end
end


%% LOG 
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
        if iRegion ==1
        violinplot([ ...
            log(d_prime{1}{iRegion}(kp0, iPair));log(d_prime{2}{iRegion}(kp1, iPair)); log(d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ones(length(d_prime{1}{iRegion}(kp0, iPair)), 1) .* 1;
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([3,1,2],:));
        else
             violinplot([ ...
          log(d_prime{2}{iRegion}(kp1, iPair)); log(d_prime{3}{iRegion}(kp2, iPair))], ...
            [ ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([1,2],:));
     
        end

        lmeTable = table;
        lmeTable.dprime = log(d_prime{2}{iRegion}(kp1, iPair));
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

        lmeTable2.dprime = log(d_prime{3}{iRegion}(kp2, iPair));
        lmeTable2.sessionType = ones(size(d_prime{3}{iRegion}(kp2, iPair), 1), 1).*2;
        lmeTable2.session = d_prime_session_num{3}{iRegion}(kp2, iPair);
        mousee = d_prime_session_num{3}{iRegion}(kp2, iPair);
        uu = unique(d_prime_session_num{3}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable2.mouse = mousee;

        if iRegion ==1
        lmeTable3 = table;

        uu = unique(d_prime_session_num{1}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{1}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable3.dprime = log(d_prime{1}{iRegion}(kp0, iPair));
        lmeTable3.sessionType = ones(size(d_prime{1}{iRegion}(kp0, iPair), 1), 1).*3;
        lmeTable3.session = d_prime_session_num{1}{iRegion}(kp0, iPair);
        mousee = d_prime_session_num{1}{iRegion}(kp0, iPair);
        uu = unique(d_prime_session_num{1}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable3.mouse = mousee;

        lmeTable13 = [lmeTable3;lmeTable2];
        lmeTable12 = [lmeTable3;lmeTable];
        end

        lmeTable23 = [lmeTable;lmeTable2];
        

        lme23 = fitlme(lmeTable23, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        if iRegion==1
        try
        lme13 = fitlme(lmeTable13, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        lme12 = fitlme(lmeTable12, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        
        catch
            lme13.coefTest=NaN;
            lme12.coefTest=NaN;
        end
        end
       
        
        %contrastMatrix = [0 1 0; 0 0 1]'; % For Type1 vs Type2, and Type1 vs Type3 comparisons
        %pValue = coefTest(lme, contrastMatrix);
        
        ylim([-10, 5])
        if iRegion==1
        prettify_pvalues(gca, [1,2,1], [2,3,3], [lme12.coefTest, lme23.coefTest, lme13.coefTest])
        else
             prettify_pvalues(gca, [1], [2], [lme23.coefTest])
       
        end


    end
end

%% x SQRT 
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
        if iRegion ==1
        violinplot([ ...
            (d_prime{1}{iRegion}(kp0, iPair)).^0.25;(d_prime{2}{iRegion}(kp1, iPair)).^0.25; (d_prime{3}{iRegion}(kp2, iPair)).^0.25], ...
            [ ones(length(d_prime{1}{iRegion}(kp0, iPair)), 1) .* 1;
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([3,1,2],:));
        else
             violinplot([ ...
          (d_prime{2}{iRegion}(kp1, iPair)).^0.25; (d_prime{3}{iRegion}(kp2, iPair)).^0.25], ...
            [ ...
            ones(length(d_prime{2}{iRegion}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{3}{iRegion}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([1,2],:));
     
        end

        lmeTable = table;
        lmeTable.dprime = (d_prime{2}{iRegion}(kp1, iPair)).^0.25;
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

        lmeTable2.dprime = (d_prime{3}{iRegion}(kp2, iPair)).^0.25;
        lmeTable2.sessionType = ones(size(d_prime{3}{iRegion}(kp2, iPair), 1), 1).*2;
        lmeTable2.session = d_prime_session_num{3}{iRegion}(kp2, iPair);
        mousee = d_prime_session_num{3}{iRegion}(kp2, iPair);
        uu = unique(d_prime_session_num{3}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable2.mouse = mousee;

        if iRegion ==1
        lmeTable3 = table;

        uu = unique(d_prime_session_num{1}{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{1}{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{1}{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable3.dprime = (d_prime{1}{iRegion}(kp0, iPair)).^0.25;
        lmeTable3.sessionType = ones(size(d_prime{1}{iRegion}(kp0, iPair), 1), 1).*3;
        lmeTable3.session = d_prime_session_num{1}{iRegion}(kp0, iPair);
        mousee = d_prime_session_num{1}{iRegion}(kp0, iPair);
        uu = unique(d_prime_session_num{1}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable3.mouse = mousee;

        lmeTable13 = [lmeTable3;lmeTable2];
        lmeTable12 = [lmeTable3;lmeTable];
        end

        lmeTable23 = [lmeTable;lmeTable2];
        

        lme23 = fitlme(lmeTable23, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        if iRegion==1
        try
        lme13 = fitlme(lmeTable13, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        lme12 = fitlme(lmeTable12, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        
        catch
            lme13.coefTest=NaN;
            lme12.coefTest=NaN;
        end
        end
       
        
        %contrastMatrix = [0 1 0; 0 0 1]'; % For Type1 vs Type2, and Type1 vs Type3 comparisons
        %pValue = coefTest(lme, contrastMatrix);
        
        ylim([-10, 5])
        if iRegion==1
        prettify_pvalues(gca, [1,2,1], [2,3,3], [lme12.coefTest, lme23.coefTest, lme13.coefTest])
        else
             prettify_pvalues(gca, [1], [2], [lme23.coefTest])
       
        end


    end
end


%% drpime of stim 1-2 throughout regions
%% SQRT 
theseColors = ya_getColors(3);
figure();
for iPair = 1%:3
    %for iRegion = 1:3
       % subplot(3, 3, iPair+(iRegion - 1)*3);
        num_recs = size(unique(d_prime_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];

        iTask = 3;
        kp0 = find(abs(d_prime{iTask}{1}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{1}(:, 1))) & ~isnan(abs(d_prime{iTask}{1}(:, 1))));
        kp1 = find(abs(d_prime{iTask}{2}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{2}(:, 1))) & ~isnan(abs(d_prime{iTask}{2}(:, 1))));
        kp2 = find(abs(d_prime{iTask}{3}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iTask}{3}(:, 1))) & ~isnan(abs(d_prime{iTask}{3}(:, 1))));


        %
        % violinplot([sqrt(d_prime{1}{iRegion}(kp0, iPair));...
        %     sqrt(d_prime{2}{iRegion}(kp1, iPair)); d_prime{3}{iRegion}(kp2, iPair)], ...
        %     [ones(length(d_prime{1}{iRegion}(kp0, iPair)),1);...
        %     ones(length(d_prime{2}{iRegion}(kp1, iPair)),1).*2; ...
        %     ones(length(d_prime{3}{iRegion}(kp2, iPair)),1).*3], ...
        %     'ViolinColor', theseColors);
       
        violinplot([ ...
            sqrt(d_prime{iTask}{1}(kp0, iPair));sqrt(d_prime{iTask}{2}(kp1, iPair)); sqrt(d_prime{iTask}{3}(kp2, iPair))], ...
            [ ones(length(d_prime{iTask}{1}(kp0, iPair)), 1) .* 1;
            ones(length(d_prime{iTask}{2}(kp1, iPair)), 1) .* 2; ...
            ones(length(d_prime{iTask}{3}(kp2, iPair)), 1) .* 3], ...
            'ViolinColor', theseColors([3,1,2],:));
      

        lmeTable = table;
        lmeTable.dprime = sqrt(d_prime{iTask}{1}(kp0, iPair));
        lmeTable.sessionType = ones(size(d_prime{iTask}{1}(kp0, iPair), 1), 1);
        lmeTable.session = d_prime_session_num{iTask}{1}(kp0, iPair);

        uu = unique(d_prime_session_num{iTask}{1});
        uu = uu(~isnan(unique(d_prime_session_num{iTask}{1}(:, iPair))));
        mouseys = [d_prime_animal_num{iTask}{1, ...
            uu(~isnan(unique(d_prime_session_num{iTask}{1}(:, iPair))), 1) ...
            }];

        mousee = d_prime_session_num{iTask}{1}(kp0, iPair);
        uu = unique(d_prime_session_num{iTask}{1});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = mousee;
        
        lmeTable2 = table;

        uu = unique(d_prime_session_num{iTask}{2});
        uu = uu(~isnan(unique(d_prime_session_num{iTask}{2}(:, iPair))));
        mouseys = [d_prime_animal_num{iTask}{2, ...
            uu(~isnan(unique(d_prime_session_num{iTask}{2}(:, iPair))), 1) ...
            }];

        lmeTable2.dprime = sqrt(d_prime{iTask}{2}(kp1, iPair));
        lmeTable2.sessionType = ones(size(d_prime{iTask}{2}(kp1, iPair), 1), 1).*2;
        lmeTable2.session = d_prime_session_num{iTask}{2}(kp1, iPair);
        mousee = d_prime_session_num{iTask}{2}(kp1, iPair);
        uu = unique(d_prime_session_num{iTask}{2});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable2.mouse = mousee;

        %if iRegion ==1
        lmeTable3 = table;

        uu = unique(d_prime_session_num{iTask}{3});
        uu = uu(~isnan(unique(d_prime_session_num{iTask}{3}(:, iPair))));
        mouseys = [d_prime_animal_num{iTask}{3, ...
            uu(~isnan(unique(d_prime_session_num{iTask}{3}(:, iPair))), 1) ...
            }];

        lmeTable3.dprime = sqrt(d_prime{iTask}{3}(kp2, iPair));
        lmeTable3.sessionType = ones(size(d_prime{iTask}{3}(kp2, iPair), 1), 1).*3;
        lmeTable3.session = d_prime_session_num{iTask}{3}(kp2, iPair);
        mousee = d_prime_session_num{iTask}{3}(kp2, iPair);
        uu = unique(d_prime_session_num{iTask}{3});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable3.mouse = mousee;

        lmeTable13 = [lmeTable;lmeTable3];
        lmeTable12 = [lmeTable;lmeTable2];
       % end

        lmeTable23 = [lmeTable2;lmeTable3];
        

        lme23 = fitlme(lmeTable23, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
       
        lme13 = fitlme(lmeTable13, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
        lme12 = fitlme(lmeTable12, 'dprime ~ sessionType + (1|mouse) + (1|mouse:session)'); %what if I add 1|neuron?
 
        
       
        
        %contrastMatrix = [0 1 0; 0 0 1]'; % For Type1 vs Type2, and Type1 vs Type3 comparisons
        %pValue = coefTest(lme, contrastMatrix);
        
        ylim([-0.5, 2.5])
        
        prettify_pvalues(gca, [1,2,1], [2,3,3], [lme12.coefTest, lme23.coefTest, lme13.coefTest])
        


    end
%end


%%  compare firing rate stim 1 and 2, pairwise 
figure();
for iTask = 1:3
%for iPair = 1:3
    for iRegion = 1:3
        subplot(3, 3, iTask+(iRegion - 1)*3);
        num_recs = size(unique(fr_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];

       % iTask = 1;
       % kp0 = find(abs(fr{iTask}{iRegion}(:, iPair)) ~= 0 & ...
       %     ~isinf(abs(fr{iTask}{iRegion}(:, 1))) & ~isnan(abs(fr{iTask}{iRegion}(:, 1))));


      %  iTask = 2;
       % kp1 = find(abs(fr{iTask}{iRegion}(:, iPair)) ~= 0 & ...
      %      ~isinf(abs(fr{iTask}{iRegion}(:, 1))) & ~isnan(abs(fr{iTask}{iRegion}(:, 1))));


        
        kp2 = find(abs(fr{iTask}{iRegion}(:, 1)) ~= 0 & ~isinf(abs(fr{iTask}{iRegion}(:, 1))) & ~isnan(abs(fr{iTask}{iRegion}(:, 1))) &...
            abs(fr{iTask}{iRegion}(:, 2)) ~= 0 & ~isinf(abs(fr{iTask}{iRegion}(:, 2))) & ~isnan(abs(fr{iTask}{iRegion}(:, 2)))&...
            abs(fr{iTask}{iRegion}(:, 3)) ~= 0 & ~isinf(abs(fr{iTask}{iRegion}(:, 3))) & ~isnan(abs(fr{iTask}{iRegion}(:, 3)))&...
            abs(fr{iTask}{iRegion}(:, 1)) < 50&abs(fr{iTask}{iRegion}(:, 2)) < 50&abs(fr{iTask}{iRegion}(:, 3)) < 50 );
hold on;
        for iCell=1:length(kp2)
            plot([1,2,3],[(fr{iTask}{iRegion}(kp2(iCell), 1));(fr{iTask}{iRegion}(kp2(iCell), 2)); (fr{iTask}{iRegion}(kp2(iCell), 3))],'Color', [0, 0, 0, 0.2])
        end
violinplot([ ...
            (fr{iTask}{iRegion}(kp2, 1));(fr{iTask}{iRegion}(kp2, 2)); (fr{iTask}{iRegion}(kp2, 3))], ...
            [ ones(length(fr{iTask}{iRegion}(kp2, 1)), 1) .* 1;
            ones(length(fr{iTask}{iRegion}(kp2, 2)), 1) .* 2; ...
            ones(length(fr{iTask}{iRegion}(kp2, 3)), 1) .* 3], ...
            'ViolinColor', theseColors([3,1,2],:));
    

        iPair =1;
        lmeTable = table;
        lmeTable.fr = [fr{iTask}{iRegion}(kp2, iPair); fr{iTask}{iRegion}(kp2, 2)];
        lmeTable.sessionType = [ones(size(fr{iTask}{iRegion}(kp2, iPair), 1), 1); ones(size(fr{iTask}{iRegion}(kp2, 2), 1), 1).*2];
        lmeTable.session = [fr_session_num{iTask}{iRegion}(kp2, iPair); fr_session_num{iTask}{iRegion}(kp2, 2)];
        lmeTable.unit = [1:size(fr_session_num{iTask}{iRegion}(kp2, iPair)), 1:size(fr_session_num{iTask}{iRegion}(kp2, iPair))]';

        uu = unique(fr_session_num{iTask}{iRegion});
        uu = uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))));
        mouseys = [fr_animal_num{iTask}{iRegion, ...
            uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))), 1) ...
            }];

        mousee = fr_session_num{iTask}{iRegion}(kp2, iPair);
        uu = unique(fr_session_num{iTask}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = [mousee; mousee];
        
        lme12 = fitlme(lmeTable, 'fr ~ sessionType + (1|mouse) + (1|mouse:session) + (1|unit)'); %what if I add 1|neuron?



         iPair =1;
        lmeTable = table;
        lmeTable.fr = [fr{iTask}{iRegion}(kp2, iPair); fr{iTask}{iRegion}(kp2, 3)];
        lmeTable.sessionType = [ones(size(fr{iTask}{iRegion}(kp2, iPair), 1), 1); ones(size(fr{iTask}{iRegion}(kp2, 3), 1), 1).*2];
        lmeTable.session = [fr_session_num{iTask}{iRegion}(kp2, iPair); fr_session_num{iTask}{iRegion}(kp2, 3)];
        lmeTable.unit = [1:size(fr_session_num{iTask}{iRegion}(kp2, iPair)), 1:size(fr_session_num{iTask}{iRegion}(kp2, iPair))]';

        uu = unique(fr_session_num{iTask}{iRegion});
        uu = uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))));
        mouseys = [fr_animal_num{iTask}{iRegion, ...
            uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))), 1) ...
            }];

        mousee = fr_session_num{iTask}{iRegion}(kp2, iPair);
        uu = unique(fr_session_num{iTask}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = [mousee; mousee];
        
        lme13 = fitlme(lmeTable, 'fr ~ sessionType + (1|mouse) + (1|mouse:session) + (1|unit)'); %what if I add 1|neuron?


         iPair =3;
        lmeTable = table;
        lmeTable.fr = [fr{iTask}{iRegion}(kp2, iPair); fr{iTask}{iRegion}(kp2, 2)];
        lmeTable.sessionType = [ones(size(fr{iTask}{iRegion}(kp2, iPair), 1), 1); ones(size(fr{iTask}{iRegion}(kp2, 2), 1), 1).*2];
        lmeTable.session = [fr_session_num{iTask}{iRegion}(kp2, iPair); fr_session_num{iTask}{iRegion}(kp2, 2)];
        lmeTable.unit = [1:size(fr_session_num{iTask}{iRegion}(kp2, iPair)), 1:size(fr_session_num{iTask}{iRegion}(kp2, iPair))]';

        uu = unique(fr_session_num{iTask}{iRegion});
        uu = uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))));
        mouseys = [fr_animal_num{iTask}{iRegion, ...
            uu(~isnan(unique(fr_session_num{iTask}{iRegion}(:, iPair))), 1) ...
            }];

        mousee = fr_session_num{iTask}{iRegion}(kp2, iPair);
        uu = unique(fr_session_num{iTask}{iRegion});
        uu =uu(~isnan(uu));
        for iS = 1:size(uu, 1)
            try
            mousee(mousee == uu(iS)) = mouseys(iS);
            catch
            end
        end
        lmeTable.mouse = [mousee; mousee];
        
        lme23 = fitlme(lmeTable, 'fr ~ sessionType + (1|mouse) + (1|mouse:session) + (1|unit)'); %what if I add 1|neuron?

        
        
       
        
        %contrastMatrix = [0 1 0; 0 0 1]'; % For Type1 vs Type2, and Type1 vs Type3 comparisons
        %pValue = coefTest(lme, contrastMatrix);
        
       % ylim([-10, 5])
        %if iRegion==1
        prettify_pvalues(gca, [1,2,1], [2,3,3], [lme12.coefTest, lme23.coefTest, lme13.coefTest])
        %else
           %  prettify_pvalues(gca, [1], [2], [lme23.coefTest])
       
        %end

    end
end
% plot
% plot histograms sorted by session and mouse (tasks next to each other)

% plot histograms sorted by session and mouse (regions next to each other)
