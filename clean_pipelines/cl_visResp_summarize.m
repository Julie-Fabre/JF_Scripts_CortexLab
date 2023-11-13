
%% cl_dprime_summarize

% parameters
keepVis = 0;
keepUnits = [1, 2]; % 1:good, 2:mua, 3:non-somatic, 4:noise
plotMe = false;
plotRegions = [1,2,3];

% initialize variables
vis_resp = cell(3, 1);
vis_resp_session_num = cell(3, 1);
vis_resp_animal_num = cell(3, 1);
vis_resp_session_fraction = cell(3, 1);
vis_resp_session_faction_median = cell(3, 1);

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

    [vis_resp{iTask}, vis_resp_session_num{iTask}, vis_resp_animal_num{iTask}, vis_resp_session_fraction{iTask}]...
        = cl_visResp(task_data, idx, keepVis, keepUnits, plotRegions, plotMe);
end

% plot overlaid histograms
cols = ya_getColors(3);
figure();
for iRegion = 1:3
    for iPair = 1:3
        for iTask = 3:-1:1
            %if iTask == 1 && (iRegion == 2 || iRegion == 3)
            %    continue;
            %end
            subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)));
            hold on;
            kp = find(abs(vis_resp{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iPair))));

            histogram(vis_resp{iTask}{iRegion}(kp, iPair), 0:0.05:4, 'Normalization', 'probability',...
                'FaceColor', cols(iTask,:), 'EdgeColor', cols(iTask,:), 'FaceAlpha', 0.1);
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
                subplot(3, 3, iPair+(iRegion- 1)*3);
                num_recs = size(unique(vis_resp_session_num{iTask}{iRegion}(:, 1)), 1);
                %theseColors = [lines(num_recs)];
                
                
                iTask = 2;
                kp1 = find(abs(vis_resp{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(vis_resp{iTask}{iRegion}(:, 1))) & ~isnan(abs(vis_resp{iTask}{iRegion}(:, 1))));
                nums1 = vis_resp_session_num{iTask}{iRegion}(kp1,1);
                [~, ~, ic] = unique(nums1, 'stable');
                nums_u1 = ic';
                
                iTask = 3;
                kp2 = find(abs(vis_resp{iTask}{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(vis_resp{iTask}{iRegion}(:, 1))) & ~isnan(abs(vis_resp{iTask}{iRegion}(:, 1))));
                nums2 = vis_resp_session_num{iTask}{iRegion}(kp2,1);
                [~, ~, ic] = unique(nums2, 'stable');
                nums_u2 = ic';

                theseColors = repelem(cols(2,:), length(unique(nums_u1)), 1);
                theseColors = [theseColors; repelem(cols(3,:), length(unique(nums_u2)), 1)];

                violinplot([vis_resp{2}{iRegion}(kp1, iPair); vis_resp{3}{iRegion}(kp2, iPair)], ...
                    [nums_u1'; max(nums_u1)+nums_u2'], ...
                    'ViolinColor', theseColors);
                ylim([0, 2.6])

                uu = unique(vis_resp_session_num{2}{iRegion});
                mouseys = [vis_resp_animal_num{2}{iRegion, ...
                    uu(~isnan(unique(vis_resp_session_num{2}{iRegion}(:, 1))), 1) ...
                    }];
                mouseys_u = unique(mouseys);
                cc =0;
                for iMousey = 1:size(mouseys_u,2)
                    idx_m = find(mouseys == mouseys_u(iMousey));
                    line([cc , cc + length(idx_m)+0.5], [4, 4], 'Color', colsMice(mouseys_u(iMousey),:));
                    cc = cc + length(idx_m);
                end
                
                uu = unique(vis_resp_session_num{3}{iRegion});
                mouseys = [vis_resp_animal_num{3}{iRegion, ...
                    uu(~isnan(unique(vis_resp_session_num{3}{iRegion}(:, 1))), 1) ...
                    }];
                mouseys_u = unique(mouseys);
                %cc =0;
                for iMousey = 1:size(mouseys_u,2)
                    idx_m = find(mouseys == mouseys_u(iMousey));
                    line([cc , cc + length(idx_m)+ 0.5], [4, 4], 'Color', colsMice(mouseys_u(iMousey),:));
                    cc = cc + length(idx_m);
                end

                axis tight; 


    end
end
prettify_plot('YLimits', 'all')

% plot histograms sorted by session and mouse (tasks next to each other)

% plot histograms sorted by session and mouse (regions next to each other)

