% dot plots, comparing firing rate in gogogo for the different stims

% (mean subtracted)

% QQ homogenize so same number conds

regions = {'CP', 'GPe', 'SNr'};
datasetlocations = {'/home/julie/Dropbox/MATLAB/naive_data4.mat', ... %6 cw stims
    '/home/julie/Dropbox/MATLAB/naive_data5.mat', ... %7 cw stims
    '/home/julie/Dropbox/MATLAB/gogogo_data2.mat', ...
    '/home/julie/Dropbox/MATLAB/goNogo_data2.mat'};
conditionsIndex = [4, 5, 2, 2];
firing_rate_per_condition = cell(3, 4);
bl_firing_rate_per_condition = cell(3, 4);
std_per_condition = cell(3, 4);
psth_per_cond = cell(3, 4);

for iDataset = 2:4

    % for each condition, measure
    passive_data_per_cond = load(datasetlocations{iDataset});

    if iDataset == 3 || iDataset == 4
        passive_data_per_cond.psth_conditions_all = passive_data_per_cond.psth_conditions_all(:, 2);
    end
    keep psth_per_cond passive_data_per_cond regions conditionsIndex iDataset datasetlocations firing_rate_per_condition bl_firing_rate_per_condition std_per_condition
    n_conditions = 1;
    %conditionsIndex = 6;%[1,2,6];
  %  figure();
    for iRegion = 1:3
        psth_all = [];

        %for iCondition = 1:n_conditions %QQ for 4th, combine 5th, + only use central ?
        [regionClassification, unitClassification] = cl_subsection_region(passive_data_per_cond.unit_coords, ...
            passive_data_per_cond.unit_area, passive_data_per_cond.pss, passive_data_per_cond.templateDuration, passive_data_per_cond.propLongISI);
        unique_subReg = unique(regionClassification(passive_data_per_cond.unit_area == iRegion));
        unique_subReg(unique_subReg == "") = [];
        %
        % if iRegion == 1
        %     curr_units = find(passive_data_per_cond.unit_area == iRegion & (passive_data_per_cond.unitType' == 1 | passive_data_per_cond.unitType' == 2) & ...
        %          unitClassification == 'MSN' & regionClassification == 'dorsomedial_striatum'& ...
        %         passive_data_per_cond.pvalue_shuffled_005{1, conditionsIndex(iDataset)}' == 1);
        % else
        %     curr_units = find(passive_data_per_cond.unit_area == iRegion & (passive_data_per_cond.unitType' == 1 | passive_data_per_cond.unitType' == 2)& ...
        %         passive_data_per_cond.pvalue_shuffled_005{1, conditionsIndex(iDataset)}' == 1);
        % end
        keepUnits = [1, 2];
        curr_units = find(passive_data_per_cond.unit_area == iRegion & ...
            ismember(passive_data_per_cond.unitType', keepUnits) & ...
            (passive_data_per_cond.pvalue_shuffled_005{1, conditionsIndex(iDataset)}(1:size(passive_data_per_cond.unit_area, 1))' == 1));

        % selectivity_index{iDataset,iRegion} = nan(size(curr_units, 1), 1);

        use_conditions = 1:size(passive_data_per_cond.psth{conditionsIndex(iDataset)}, 3);

        nonEmptyRecs = find(~cellfun(@isempty, passive_data_per_cond.psth_conditions_all));
        passive_data_per_cond.animal_thisDate_site_shank(isnan(passive_data_per_cond.animal_thisDate_site_shank(:, 4)), 4) = 0;

        unique_recs = unique(passive_data_per_cond.animal_thisDate_site_shank, 'rows');

        for iUnit = 1:size(curr_units, 1)

            thisUnit = curr_units(iUnit);
            [~, iRec] = ismember(passive_data_per_cond.animal_thisDate_site_shank(thisUnit, :), unique_recs, 'rows');
            try
                conditions = passive_data_per_cond.psth_conditions_all{nonEmptyRecs(iRec)};
            catch
                continue;
            end
            if sum(nanmean(squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit, 1, :, 250:450))) == 0) > 190
                continue;
            end

            if iDataset == 1 % cw stims, stim 1 v 2
                %conditions = conditions(:,1);
                use_conditions = conditions(ismember(conditions(:, 1), [1, 2, 3]) & conditions(:, 2) == -90, :);
                [~, condType] = ismember(conditions, use_conditions, 'rows');
            elseif iDataset == 2 % cw stims, stim 1 v 2
                %conditions = conditions(:,1);
                if sum(conditions(:, 1) == 12) > 0
                    use_conditions = conditions(ismember(conditions(:, 1), [4, 6, 12]) & conditions(:, 2) == -90, :);
                else
                    use_conditions = conditions(ismember(conditions(:, 1), [1, 2, 3]) & conditions(:, 2) == -90, :);
                end
                [~, condType] = ismember(conditions, use_conditions, 'rows');
            elseif iDataset == 3 || iDataset == 4 % cw stims in tasks, stim 1 v 2
                if sum(conditions(:, 1) == 12) > 0
                    try
                        use_conditions = conditions(ismember(conditions(:, 1), [4, 6, 12]) & conditions(:, 2) == -90, :);
                    catch
                        keyboard
                    end
                else
                    try
                        use_conditions = conditions(ismember(conditions(:, 1), [1, 2, 3]) & conditions(:, 2) == -90, :);
                    catch
                        keyboard
                    end
                end
                [~, condType] = ismember(conditions, use_conditions, 'rows');

            end


            %condIdx = [17, 19, 25];
            if size(find(condType)) ~= 3
                continue;
            end
            for iCond = 1:size(use_conditions, 1)
                baseline_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit, ...
                    3, condType == iCond, 1:100), 3), 4)) .* 100;
                resp_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit, ...
                    3, condType == iCond, 260:360), 3), 4)) .* 100;
                std_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit, ...
                    3, condType == iCond, 1:100), 3), 4)) .* 100;
                psth_all(iUnit, iCond, :) = squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit, ...
                    3, condType == iCond, :));
            end

            % % DEBUGGING
            % figure();
            % plot(smoothdata(squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
            %         3, condType == 1, :)), 'gaussian', [0 50]))
            % hold on;
            % plot(smoothdata(squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
            %         3, condType == 2, :)), 'gaussian', [0 50]))
            % plot(smoothdata(squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
            %         3, condType == 3, :)), 'gaussian', [0 50]))
            % hold off;
            % title(num2str(resp_per_cond))

            firing_rate_per_condition{iRegion, iDataset}(iUnit, :) = resp_per_cond;
            bl_firing_rate_per_condition{iRegion, iDataset}(iUnit, :) = baseline_per_cond;
            std_per_condition{iRegion, iDataset}(iUnit, :) = std_per_cond;
        end

        % subplot(3, 3, (3 * (iRegion - 1))+1)
        % %subplot(131)
        smooth_filt = [1, 10]; % (units x frames)
        img1 = (squeeze(psth_all(:, 1, :)) - nanmean(psth_all(:, 1, 1:100), 3)) ./ (nanstd(psth_all(:, 1, 1:100), [], 3));
        [~, cell_idx] = sort(nanmean(img1(:, 260:360), 2));
        img1_sm = conv2(img1(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img1(cell_idx, :)), ones(smooth_filt), 'same');
        % imagesc(img1_sm)
        % caxis([-4 4])
        %
        % subplot(3, 3, (3 * (iRegion - 1))+3)
        img2 = (squeeze(psth_all(:, 2, :)) - nanmean(psth_all(:, 2, 1:100), 3)) ./ (nanstd(psth_all(:, 2, 1:100), [], 3));
        img2_sm = conv2(img2(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img2(cell_idx, :)), ones(smooth_filt), 'same');
        % imagesc(img2_sm)
        % caxis([-4 4])
        %
        % subplot(3, 3, (3 * (iRegion - 1))+2)
        img3 = (squeeze(psth_all(:, 3, :)) - nanmean(psth_all(:, 3, 1:100), 3)) ./ (nanstd(psth_all(:, 3, 1:100), [], 3));
        img3_sm = conv2(img3(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img3(cell_idx, :)), ones(smooth_filt), 'same');
        % imagesc(img3_sm)
        % caxis([-4 4])

        colormap(brewermap([], '*RdBu'))


        psth_per_cond{iRegion, iDataset}(:, :, 1) = img1_sm;
        psth_per_cond{iRegion, iDataset}(:, :, 2) = img2_sm;
        psth_per_cond{iRegion, iDataset}(:, :, 3) = img3_sm;

    end
end

%% correlation
dot_size = 3;
for iDataset = 3:4
    figure();
    for iRegion = 1:3
        % remove zero rows
        nonEmpty_rows = sum(firing_rate_per_condition{iRegion, iDataset}, 2) > 0;
        % plot

        act_go2 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 3); %go2
        act_stim3 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 2); %nogo
        act_go1 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 1); %go1

        subplot(3, 1, iRegion)
        % Determine the number of neurons
        numNeurons = size(psth_per_cond{iRegion, iDataset}, 1);

        % Initialize a vector to hold the correlation coefficients
        correlationCoefficients = zeros(3, numNeurons, 1);

        for iNeuron = 1:numNeurons
            % Extract the activity for the current neuron in both conditions
            activityCond1 = psth_per_cond{iRegion, iDataset}(iNeuron, :, 3);
            activityCond2 = psth_per_cond{iRegion, iDataset}(iNeuron, :, 1);
            activityCond3 = psth_per_cond{iRegion, iDataset}(iNeuron, :, 2);


            % Calculate the Pearson correlation coefficient for this neuron
            % between the two conditions. Note: corr requires column vectors,
            % so ensure the data is appropriately shaped.
            R12 = corr(activityCond1(:), activityCond2(:));
            R23 = corr(activityCond2(:), activityCond3(:));
            R13 = corr(activityCond1(:), activityCond3(:));

            % Store the coefficient
            correlationCoefficients(1, iNeuron) = R12;
            correlationCoefficients(2, iNeuron) = R23;
            correlationCoefficients(3, iNeuron) = R13;
        end

        boxplot([correlationCoefficients(1, :), correlationCoefficients(2, :), correlationCoefficients(3, :)], ...
            [ones(numNeurons, 1); ones(numNeurons, 1) .* 2; ones(numNeurons, 1) .* 3])
        xticks([1, 2, 3])
        xticklabels({'go1 v go2', 'go2 v stim 3', 'go1 v stim 3'})
        ylabel('correlation coefficient')


    end
    prettify_plot('YLimits', 'all', 'XLimits', 'all');
end

%% firing rate dot plot
dot_size = 3;
texty = [50,4,1];
for iDataset = 3:4
    figure();
    for iRegion = 1:3
        % remove zero rows
        nonEmpty_rows = sum(firing_rate_per_condition{iRegion, iDataset}, 2) > 0;
        % plot
        subplot(3, 3, (iRegion - 1)*3+1)
        act_go2 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 3); %go2
        act_stim3 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 2); %nogo
        act_go1 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 1); %go1

        indexes = ~isinf(act_go2) & ~isinf(act_go1) & ~isinf(act_stim3) & ~isnan(act_go2) & ~isnan(act_go1) & ~isnan(act_stim3);


        scatter(act_go2(indexes), act_go1(indexes), dot_size, 'filled');
        hold on;
        %[rho, pval] = corr(act_go1, act_go2);
        p = polyfit(act_go2(indexes), act_go1(indexes), 1); %1 degree: linear
        slope = p(1);
        x_fit = linspace(min(act_go2(indexes)), max(act_go2(indexes)), 100); % Generate 100 points between the min and max of act_go1
        y_fit = polyval(p, x_fit); % Evaluate the polynomial p at points in x_fit
        % Plot the linear regression line
        %plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % 'r-' for a red line, 'LineWidth' to define the thickness
        %x_text = texty(iRegion); %max(x_fit); % x position for the text
        %y_text = polyval(p, x_text); % y position for the text, calculated using the regression line equation
        %textStr = ['Slope: ', num2str(p(1), '%.2f')]; % Text string with slope value rounded to 2 decimal places
        %text(x_text, y_text, textStr, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        xlabel('go 2 activity (zscored)')
        ylabel('go 1 activity (zscored)')
         line([-1, texty(iRegion)],[-1, texty(iRegion)])
        title(regions{iRegion})
        hold off;

        % xlim([0.001, 100])
        % ylim([0.001, 100])
        % line([0.001, 100],[0.001, 100])
        % set(gca, 'XScale', 'log', 'YScale', 'log');

        subplot(3, 3, (iRegion - 1)*3+2)
        scatter(act_go2(indexes), act_stim3(indexes), dot_size, 'filled');
        hold on;
        p = polyfit(act_go2(indexes), act_stim3(indexes), 1); %1 degree: linear
        slope = p(1);
        x_fit = linspace(min(act_go2(indexes)), max(act_go2(indexes)), 100); % Generate 100 points between the min and max of act_go1
        y_fit = polyval(p, x_fit); % Evaluate the polynomial p at points in x_fit
        % Plot the linear regression line
        %plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % 'r-' for a red line, 'LineWidth' to define the thickness
        %x_text = texty(iRegion); % x position for the text
        %y_text = polyval(p, x_text); % y position for the text, calculated using the regression line equation
        %textStr = ['Slope: ', num2str(p(1), '%.2f')]; % Text string with slope value rounded to 2 decimal places
        %text(x_text, y_text, textStr, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        xlabel('go 2 activity (zscored)')
        ylabel('stim 3 activity (zscored)')
         line([-1, texty(iRegion)],[-1, texty(iRegion)])
        title(regions{iRegion})
        hold off;

        subplot(3, 3, (iRegion - 1)*3+3)
        scatter(act_go1(indexes), act_stim3(indexes), dot_size, 'filled');
        hold on;
        p = polyfit(act_go1(indexes), act_stim3(indexes), 1); %1 degree: linear
        slope = p(1);
        x_fit = linspace(min(act_go1(indexes)), max(act_go1(indexes)), 100); % Generate 100 points between the min and max of act_go1
        y_fit = polyval(p, x_fit); % Evaluate the polynomial p at points in x_fit
        % Plot the linear regression line
        %plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % 'r-' for a red line, 'LineWidth' to define the thickness
        %x_text = texty(iRegion); % x position for the text
        %y_text = polyval(p, x_text); % y position for the text, calculated using the regression line equation
        %textStr = ['Slope: ', num2str(p(1), '%.2f')]; % Text string with slope value rounded to 2 decimal places
        %text(x_text, y_text, textStr, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        hold off;
        % xlim([0.001, 100])
        % ylim([0.001, 100])
        line([-1, texty(iRegion)],[-1,texty(iRegion)])
        % set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('go 1 activity (zscored)')
        ylabel('stim 3 activity (zscored)')
        title(regions{iRegion})
        hold off;
    end

    prettify_plot('YLimits', 'rows', 'XLimits', 'rows')
end

prettify_plot('YLimits', [-1, texty(1)], 'XLimits', [-1, texty(1)])
prettify_plot('YLimits', [-1, texty(2)], 'XLimits', [-1, texty(2)])
prettify_plot('YLimits', [-1, texty(3)], 'XLimits', [-1, texty(3)])
%% paired firing rate plot

%% naive vs gogogo/gonogo: resp. to stimuli

for iDataset_comp = 3:4
    figure();
    for iRegion = 1:3
    
        
    iDataset = 2;
        nonEmpty_rows = sum(firing_rate_per_condition{iRegion, iDataset}, 2) > 0;

        psth_per_cond{iRegion, iDataset}(:, :, 1)
        act_go2_n = abs((firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)); %go2
        act_stim3_n = abs((firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)); %nogo
        act_go1_n =abs( (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)); %go1

        indexes_n = ~isinf(act_go2_n) & ~isinf(act_go1_n) & ~isinf(act_stim3_n) & ~isnan(act_go2_n) & ~isnan(act_go1_n) & ~isnan(act_stim3_n);

         iDataset = iDataset_comp;
        nonEmpty_rows = sum(firing_rate_per_condition{iRegion, iDataset}, 2) > 0;

        psth_per_cond{iRegion, iDataset}(:, :, 1)
        act_go2 = abs((firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 3)); %go2
        act_stim3 = abs((firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 2)); %nogo
        act_go1 = abs((firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1) - bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)) ./ ...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows, 1)); %go1

        indexes = ~isinf(act_go2) & ~isinf(act_go1) & ~isinf(act_stim3) & ~isnan(act_go2) & ~isnan(act_go1) & ~isnan(act_stim3);


        subplot(3, 3, (iRegion - 1)*3+1)
        boxplot([act_go1_n(indexes_n); act_go1(indexes)], ...
            [ones(length(act_go1_n(indexes_n)), 1); ones(length(act_go1(indexes)), 1) .* 2])
        xticks([1, 2])
        if iDataset_comp==3
            xticklabels({'naive', 'gogogo'})
        else
            xticklabels({'naive', 'go-Nogo'})
        end

        subplot(3, 3, (iRegion - 1)*3+2)
        boxplot([act_go2_n(indexes_n); act_go2(indexes)], ...
            [ones(length(act_go2_n(indexes_n)), 1); ones(length(act_go2(indexes)), 1) .* 2])
        xticks([1, 2])
        if iDataset_comp==3
            xticklabels({'naive', 'gogogo'})
        else
            xticklabels({'naive', 'go-Nogo'})
        end

        subplot(3, 3, (iRegion - 1)*3+3)
        boxplot([act_stim3_n(indexes_n); act_stim3(indexes)], ...
            [ones(length(act_stim3_n(indexes_n)), 1); ones(length(act_stim3(indexes)), 1) .* 2])
        xticks([1, 2])
        if iDataset_comp==3
            xticklabels({'naive', 'gogogo'})
        else
            xticklabels({'naive', 'go-Nogo'})
        end
    end
end
    

prettify_plot('YLimits', 'rows', 'XLimits', 'rows')


for iP=1:9
    subplot(3, 3, iP)
    %xlim([-5 50])
    ylim([-0.1 5])
end
%% naive vs goNogo: resp. to stimuli