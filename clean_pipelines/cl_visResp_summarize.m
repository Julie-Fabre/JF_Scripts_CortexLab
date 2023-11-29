
%% cl_dprime_summarize

% parameters
keepVis = 0;
keepUnits = [1]; %, 2]; % 1:good, 2:mua, 3:non-somatic, 4:noise
plotMe = false;
plotRegions = [1, 2, 3];

% initialize variables
vis_resp = cell(3, 1);
vis_resp_session_num = cell(3, 1);
vis_resp_animal_num = cell(3, 1);
vis_resp_session_fraction = cell(3, 1);
vis_resp_session_faction_median = cell(3, 1);

% run loop
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full2.mat');
        idx = 5;
        passive = 1;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo2.mat');
        idx = 2;
        passive = 0;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo2.mat');
        idx = 2;
        passive = 0;

    end
    [vis_resp{iTask}, vis_resp1{iTask}, vis_resp2{iTask}, vis_resp_session_num{iTask}, vis_resp_animal_num{iTask}, vis_resp_session_fraction{iTask}, ...
        vis_resp_full{iTask}, vis_resp_full_1{iTask}, vis_resp_full_2{iTask}, pvalue_shuff{iTask},pvalue_glmTest{iTask}, zscore_val{iTask}, zscore_tiny_tiny{iTask}, ...
        zscore_tiny{iTask}, dfr{iTask}, dsum_diff{iTask}] = cl_visResp(task_data, idx, keepVis, keepUnits, plotRegions, plotMe);
end

%for iTask = 1:3
    figure();
    for iPair = 1:3
        for iRegion = 1:3
            subplot(3, size(plotRegions, 2), iPair+(iRegion - 1)*(size(plotRegions, 2)));
            hold on;

            violinplot([squeeze(vis_resp_session_fraction{1}(iRegion,:,iPair)), squeeze(vis_resp_session_fraction{2}(iRegion,:,iPair)),...
                squeeze(vis_resp_session_fraction{3}(iRegion,:,iPair))], ...
                [ones(length(squeeze(vis_resp_session_fraction{1}(iRegion,:,iPair))),1);...
                ones(length(squeeze(vis_resp_session_fraction{2}(iRegion,:,iPair))),1).*2;...
                ones(length(squeeze(vis_resp_session_fraction{3}(iRegion,:,iPair))),1).*3])%iR, iS, iP
            p12 = ranksum(squeeze(vis_resp_session_fraction{1}(iRegion,:,iPair)), squeeze(vis_resp_session_fraction{2}(iRegion,:,iPair)));
            p13 = ranksum(squeeze(vis_resp_session_fraction{1}(iRegion,:,iPair)), squeeze(vis_resp_session_fraction{3}(iRegion,:,iPair)));
            p23 = ranksum(squeeze(vis_resp_session_fraction{2}(iRegion,:,iPair)), squeeze(vis_resp_session_fraction{3}(iRegion,:,iPair)));
            ylim([-0.01, 0.1])
            prettify_pvalues(gca, [1,1,2], [2,3,3], [p12, p13, p23])
            
        end

    end
%end
% sorted by increase. cross validated.
for iTask = 1:3
    figure();
    for iPair = 1:3
        for iRegion = 1:3
            %if iTask == 1 && (iRegion == 2 || iRegion == 3)
            %    continue;
            %end
            subplot(3, size(plotRegions, 2), iPair+(iRegion - 1)*(size(plotRegions, 2)));
            hold on;

            %vis resp
            kp = ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 1))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 1))) & ...
                ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 2))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 2))) & ...
                ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 3))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 3))) & ...
                (pvalue_shuff{iTask}{iRegion}(:, 1) == 1 | pvalue_shuff{iTask}{iRegion}(:, 1) ==0  | ...
               pvalue_shuff{iTask}{iRegion}(:, 2) == 1 | pvalue_shuff{iTask}{iRegion}(:, 2) ==0  | ...
                pvalue_shuff{iTask}{iRegion}(:, 3) == 1  | pvalue_shuff{iTask}{iRegion}(:, 3) ==0 ) &...
                sum(vis_resp_full_1{iTask}{iRegion}(:, 1, :),3) ~= 0 &...
                sum(vis_resp_full_1{iTask}{iRegion}(:, 2, :),3) ~= 0 &...
                sum(vis_resp_full_1{iTask}{iRegion}(:, 3, :),3) ~= 0 &...
                ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 1, :),3)) &...
                ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 2, :),3)) &...
                ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 3, :),3));


            % %pval shuff 
            %  kp = ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 1))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 1))) & ...
            %     ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 2))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 2))) & ...
            %     ~isnan(abs(vis_resp1{iTask}{iRegion}(:, 3))) & ~isinf(abs(vis_resp1{iTask}{iRegion}(:, 3))) & ...
            %     (pvalue_shuff{iTask}{iRegion}(:, 1) == 1  | ...
            %     pvalue_shuff{iTask}{iRegion}(:, 2) == 1  | ...
            %     pvalue_shuff{iTask}{iRegion}(:, 3) == 1)  &...
            %     sum(vis_resp_full_1{iTask}{iRegion}(:, 1, :),3) ~= 0 &...
            %     sum(vis_resp_full_1{iTask}{iRegion}(:, 2, :),3) ~= 0 &...
            %     sum(vis_resp_full_1{iTask}{iRegion}(:, 3, :),3) ~= 0 &...
            %     ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 1, :),3)) &...
            %     ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 2, :),3)) &...
            %     ~isnan(sum(vis_resp_full_1{iTask}{iRegion}(:, 3, :),3));


            [~, idx] = sort(vis_resp2{iTask}{iRegion}(kp,1));
            ss = squeeze(vis_resp_full_1{iTask}{iRegion}(kp, iPair, :));
            this_image = (ss(idx, :) - nanmean(ss(idx, 1:200), 2)) ./ (nanmean(ss(idx, 1:200), 2) +0.0001);
            %remove any zero rows 
            smooth_filt = [round(size(ss, 1)/20), 50];
            this_image_smooth = conv2(this_image, ones(smooth_filt), 'same') ./ ...
                conv2(~isnan(this_image), ones(smooth_filt), 'same');

            imagesc(this_image_smooth)
            %xlabel('(resp - resp0)/(resp + resp0)')
            %xlabel('(resp - resp0)/(std_resp + 0.001)')
            if iTask == 1 && iPair == 1 && iRegion == 1
                ylabel('neurons, sorted by vis')


            end
            colorbar;
            if iRegion == 1
                clim([-1, 1])
            elseif iRegion == 2
                clim([-1, 1])
            elseif iRegion == 3
                clim([-1, 1])
            end
            colormap(brewermap([], '*RdBu'))
        end
    end
    prettify_plot('AxisAspectRatio', 'tight');
end


% % sorted by the value we are interested in (not cross-val)
% %vis_resp = dfr;
% for iTask = 1:3
%     figure();
%     for iPair = 1:3
%         for iRegion = 1:3
%             %if iTask == 1 && (iRegion == 2 || iRegion == 3)
%             %    continue;
%             %end
%             subplot(3, size(plotRegions, 2), iPair+(iRegion - 1)*(size(plotRegions, 2)));
%             hold on;
%             kp = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iPair))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iPair))) & ...
%                 (vis_resp{iTask}{iRegion}(:, 1) >= 0.15 | vis_resp{iTask}{iRegion}(:, 1) <= -0.15 | ...
%                 vis_resp{iTask}{iRegion}(:, 2) >= 0.15 | vis_resp{iTask}{iRegion}(:, 2) <= -0.15 | ...
%                 vis_resp{iTask}{iRegion}(:, 3) >= 0.15 | vis_resp{iTask}{iRegion}(:, 3) <= -0.15));
% 
%             [~, idx] = sort(vis_resp{iTask}{iRegion}(kp, 2));
%             ss = squeeze(vis_resp_full{iTask}{iRegion}(kp, iPair, :));
%             imagesc(smoothdata((ss(idx, :) - nanmean(ss(idx, 1:200), 2))./(nanmean(ss(idx, 1:200), 2)), 2, 'gaussian', [50, 0]))
%             %xlabel('(resp - resp0)/(resp + resp0)')
%             %xlabel('(resp - resp0)/(std_resp + 0.001)')
%             ylabel('neurons, sorted by vis')
%             clim([-1.5, 1.5])
%             colormap(brewermap([], '*RdBu'))
%             colorbar;
% 
%         end
%     end
%     prettify_plot('AxisAspectRatio', 'tight');
% end

vis_resp =  dsum_diff
% plot overlaid histograms
cols = ya_getColors(3);
figure();
%vals = [-1:0.05:1];
vals = [-2:0.1:2];
for iRegion = 1:3
    for iStim = 1:3
        for iTask = 3:-1:1
            %if iTask == 1 && (iRegion == 2 || iRegion == 3)
            %    continue;
            %end
            subplot(3, size(plotRegions, 2), iStim+(iRegion - 1)*(size(plotRegions, 2)));
            hold on;
            kp = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))) &...
                pvalue_shuff{iTask}{iRegion}(:, 1) == 1 | pvalue_shuff{iTask}{iRegion}(:, 2) == 1 | pvalue_shuff{iTask}{iRegion}(:, 3) == 1);


            [counts, edges] = histcounts((vis_resp{iTask}{iRegion}(kp, iStim)), vals, 'Normalization', 'probability');
            stairs(edges(1:end-1), counts, 'Color', cols(iTask, :), 'LineWidth', 1.5);


            %xlabel('(resp - resp0)/(resp + resp0)')
            xlabel('(resp - resp0)/(resp + resp0)')
            ylabel('frac. neurons')
            %xlim([0, 2])
            if iTask == 1 && iRegion == 1 && iStim == 1
                legend({'go-nogo', 'gogogo', 'naive'})
            end
        end
    end
end
prettify_plot('YLimits', 'all')
 % bar plot with threshold 
figure()
for iStim = 1:3
    for iRegion = 1:3
        subplot(3, size(plotRegions, 2), iStim+(iRegion - 1)*(size(plotRegions, 2)));

        iTask = 1;

        kp1 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 2;

        kp2 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 3;

        kp3 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));

        violinplot([abs(vis_resp{1}{iRegion}(kp1, iStim))>=0.04; abs(vis_resp{2}{iRegion}(kp2, iStim))>=0.04; ...
            abs(vis_resp{3}{iRegion}(kp3, iStim))>=0.04], ...
            [ones(length(vis_resp{1}{iRegion}(kp1, iStim)), 1); ones(length(vis_resp{2}{iRegion}(kp2, iStim)), 1) .* 2; ...
            ones(length(vis_resp{3}{iRegion}(kp3, iStim)), 1) .* 3]);
        % [p1, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{2}{iRegion}(kp2, iStim));
        % [p2, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        % [p3, ~] = ranksum(vis_resp{2}{iRegion}(kp2, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        % prettify_pvalues(gca, [1,1,2], [2,3,3], [p1, p2, p3])

    end
end
prettify_plot('YLimits', 'all')

% stat test?
figure()
for iStim = 1:3
    for iRegion = 1:3
        subplot(3, size(plotRegions, 2), iStim+(iRegion - 1)*(size(plotRegions, 2)));

        iTask = 1;

        kp1 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 2;

        kp2 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 3;

        kp3 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));

        violinplot([abs(vis_resp{1}{iRegion}(kp1, iStim)); abs(vis_resp{2}{iRegion}(kp2, iStim)); ...
            abs(vis_resp{3}{iRegion}(kp3, iStim))], ...
            [ones(length(vis_resp{1}{iRegion}(kp1, iStim)), 1); ones(length(vis_resp{2}{iRegion}(kp2, iStim)), 1) .* 2; ...
            ones(length(vis_resp{3}{iRegion}(kp3, iStim)), 1) .* 3]);
        % [p1, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{2}{iRegion}(kp2, iStim));
        % [p2, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        % [p3, ~] = ranksum(vis_resp{2}{iRegion}(kp2, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        % prettify_pvalues(gca, [1,1,2], [2,3,3], [p1, p2, p3])

    end
end
prettify_plot('YLimits', 'all')

for iStim = 1:3
    for iRegion = 1:3
        subplot(3, size(plotRegions, 2), iStim+(iRegion - 1)*(size(plotRegions, 2)));
        iTask = 1;

        kp1 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 2;

        kp2 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        iTask = 3;

        kp3 = find(~isnan(abs(vis_resp{iTask}{iRegion}(:, iStim))) & ~isinf(abs(vis_resp{iTask}{iRegion}(:, iStim))));
        lmeTable = table;
        lmeTable.value = [abs(vis_resp{1}{iRegion}(kp1, iStim)); abs(vis_resp{2}{iRegion}(kp2, iStim))];
        lmeTable.group = [ones(length(vis_resp{1}{iRegion}(kp1, iStim)), 1); ones(length(vis_resp{2}{iRegion}(kp2, iStim)), 1) .* 2];
        lmeTable.session = [vis_resp_session_num{1}{iRegion}(kp1, iStim); vis_resp_session_num{2}{iRegion}(kp2, iStim)];
        %lmeTable.mouse = [vis_resp_animal_num{1}{iRegion}(kp1, iStim); vis_resp_animal_num{2}{iRegion}(kp2, iStim)];
        ff = fitlme(lmeTable, 'value ~ group + (1|session)');

        lmeTable = table;
        lmeTable.value = [abs(vis_resp{1}{iRegion}(kp1, iStim)); abs(vis_resp{3}{iRegion}(kp3, iStim))];
        lmeTable.group = [ones(length(vis_resp{1}{iRegion}(kp1, iStim)), 1); ones(length(vis_resp{3}{iRegion}(kp3, iStim)), 1) .* 2];
        lmeTable.session = [vis_resp_session_num{1}{iRegion}(kp1, iStim); vis_resp_session_num{3}{iRegion}(kp3, iStim)];
        %lmeTable.mouse = [vis_resp_animal_num{1}{iRegion}(kp1, iStim); vis_resp_animal_num{2}{iRegion}(kp2, iStim)];
        ff13 = fitlme(lmeTable, 'value ~ group + (1|session)');

        lmeTable = table;
        lmeTable.value = [abs(vis_resp{2}{iRegion}(kp2, iStim)); abs(vis_resp{3}{iRegion}(kp3, iStim))];
        lmeTable.group = [ones(length(vis_resp{2}{iRegion}(kp2, iStim)), 1); ones(length(vis_resp{3}{iRegion}(kp3, iStim)), 1) .* 2];
        lmeTable.session = [vis_resp_session_num{2}{iRegion}(kp2, iStim); vis_resp_session_num{3}{iRegion}(kp3, iStim)];
        %lmeTable.mouse = [vis_resp_animal_num{1}{iRegion}(kp1, iStim); vis_resp_animal_num{2}{iRegion}(kp2, iStim)];
        ff23 = fitlme(lmeTable, 'value ~ group + (1|session)');

        %ff.coefTest

        %[p1, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{2}{iRegion}(kp2, iStim));
        %[p2, ~] = ranksum(vis_resp{1}{iRegion}(kp1, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        %[p3, ~] = ranksum(vis_resp{2}{iRegion}(kp2, iStim), vis_resp{3}{iRegion}(kp3, iStim));
        prettify_pvalues(gca, [1, 1, 2], [2, 3, 3], [ff.coefTest, ff13.coefTest, ff23.coefTest])

    end
end
% plot hi

% plot histograms sorted by session and mouse (pairs next to each other)
figure();
for iStim = 1:3
    for iRegion = 1:3
        subplot(3, 3, iStim+(iRegion - 1)*3);
        num_recs = size(unique(vis_resp_session_num{iTask}{iRegion}(:, 1)), 1);
        %theseColors = [lines(num_recs)];


        iTask = 2;
        kp1 = find(abs(vis_resp{iTask}{iRegion}(:, iStim)) ~= 0 & ~isinf(abs(vis_resp{iTask}{iRegion}(:, 1))) & ~isnan(abs(vis_resp{iTask}{iRegion}(:, 1))));
        nums1 = vis_resp_session_num{iTask}{iRegion}(kp1, 1);
        [~, ~, ic] = unique(nums1, 'stable');
        nums_u1 = ic';

        iTask = 3;
        kp2 = find(abs(vis_resp{iTask}{iRegion}(:, iStim)) ~= 0 & ~isinf(abs(vis_resp{iTask}{iRegion}(:, 1))) & ~isnan(abs(vis_resp{iTask}{iRegion}(:, 1))));
        nums2 = vis_resp_session_num{iTask}{iRegion}(kp2, 1);
        [~, ~, ic] = unique(nums2, 'stable');
        nums_u2 = ic';

        theseColors = repelem(cols(2, :), length(unique(nums_u1)), 1);
        theseColors = [theseColors; repelem(cols(3, :), length(unique(nums_u2)), 1)];

        violinplot([vis_resp{2}{iRegion}(kp1, iStim); vis_resp{3}{iRegion}(kp2, iStim)], ...
            [nums_u1'; max(nums_u1) + nums_u2'], ...
            'ViolinColor', theseColors);
        ylim([0, 2.6])

        uu = unique(vis_resp_session_num{2}{iRegion});
        mouseys = [vis_resp_animal_num{2}{iRegion, ...
            uu(~isnan(unique(vis_resp_session_num{2}{iRegion}(:, 1))), 1)}];
        mouseys_u = unique(mouseys);
        cc = 0;
        colsMice = ya_getColors(max(mouseys_u));
        for iMousey = 1:size(mouseys_u, 2)
            idx_m = find(mouseys == mouseys_u(iMousey));
            line([cc, cc + length(idx_m) + 0.5], [4, 4], 'Color', colsMice(mouseys_u(iMousey), :));
            cc = cc + length(idx_m);
        end

        uu = unique(vis_resp_session_num{3}{iRegion});
        mouseys = [vis_resp_animal_num{3}{iRegion, ...
            uu(~isnan(unique(vis_resp_session_num{3}{iRegion}(:, 1))), 1)}];
        mouseys_u = unique(mouseys);
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

% plot histograms sorted by session and mouse (tasks next to each other)

% plot histograms sorted by session and mouse (regions next to each other)
