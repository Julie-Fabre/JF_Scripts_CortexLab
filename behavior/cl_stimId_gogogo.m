% dot plots, comparing firing rate in gogogo for the different stims 

% (mean subtracted)

% QQ homogenize so same number conds 

regions = {'CP', 'GPe', 'SNr'};
datasetlocations = { '/home/julie/Dropbox/MATLAB/naive_data4.mat', ...%6 cw stims
    '/home/julie/Dropbox/MATLAB/naive_data5.mat', ...%7 cw stims
    '/home/julie/Dropbox/MATLAB/gogogo_data2.mat', ...
    '/home/julie/Dropbox/MATLAB/goNogo_data2.mat'};
conditionsIndex = [ 4, 5, 2, 2];
firing_rate_per_condition = cell(3,4);
bl_firing_rate_per_condition = cell(3,4);
std_per_condition = cell(3,4);
psth_per_cond = cell(3,4);

for iDataset = 3:4

    % for each condition, measure
    passive_data_per_cond = load(datasetlocations{iDataset});
    
    if iDataset == 3 || iDataset == 4
            passive_data_per_cond.psth_conditions_all = passive_data_per_cond.psth_conditions_all(:,2);
    end
    keep psth_per_cond passive_data_per_cond regions conditionsIndex iDataset datasetlocations  firing_rate_per_condition bl_firing_rate_per_condition std_per_condition
    n_conditions = 1;
    %conditionsIndex = 6;%[1,2,6];
figure();
    for iRegion = 1:3
    psth_all =[];
    
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
        keepUnits = [1,2];
        curr_units = find(passive_data_per_cond.unit_area == iRegion & ...
          ismember(passive_data_per_cond.unitType', keepUnits) &...
    (passive_data_per_cond.pvalue_shuffled_005{1,conditionsIndex(iDataset)}(1:size(passive_data_per_cond.unit_area, 1))' == 1));

       % selectivity_index{iDataset,iRegion} = nan(size(curr_units, 1), 1);

        use_conditions = 1:size(passive_data_per_cond.psth{conditionsIndex(iDataset)}, 3);
        
        nonEmptyRecs = find(~cellfun(@isempty, passive_data_per_cond.psth_conditions_all));
        passive_data_per_cond.animal_thisDate_site_shank(isnan(passive_data_per_cond.animal_thisDate_site_shank(:,4)),4)=0;
       
        unique_recs = unique(passive_data_per_cond.animal_thisDate_site_shank, 'rows');

        for iUnit = 1:size(curr_units, 1)

            thisUnit = curr_units(iUnit);
            [~, iRec]= ismember(passive_data_per_cond.animal_thisDate_site_shank(thisUnit,:), unique_recs, 'rows');
            try
                conditions = passive_data_per_cond.psth_conditions_all{nonEmptyRecs(iRec)};
            catch
                continue;
            end
            if sum(nanmean(squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,1,:,250:450)))==0)>190
                continue;
            end

            if iDataset == 1 % cw stims, stim 1 v 2
                %conditions = conditions(:,1);
                use_conditions = conditions(ismember(conditions(:,1), [1,2,3]) & conditions(:,2)==-90,:);
                [~, condType] = ismember(conditions, use_conditions, 'rows');
            elseif iDataset == 2 % cw stims, stim 1 v 2
                %conditions = conditions(:,1);
                if sum(conditions(:,1) == 12) > 0
                    use_conditions = conditions(ismember(conditions(:,1), [4,6,12]) & conditions(:,2)==-90,:);
                else
                    use_conditions = conditions(ismember(conditions(:,1), [1,2,3]) & conditions(:,2)==-90,:);
                end
                [~, condType] = ismember(conditions, use_conditions, 'rows');
            elseif iDataset == 3 || iDataset == 4 % cw stims in tasks, stim 1 v 2
                if sum(conditions(:,1) == 12) > 0
                    try
                    use_conditions = conditions(ismember(conditions(:,1), [4,6,12]) & conditions(:,2)==-90,:);
                    catch
                        keyboard
                    end
                else
                    try
                    use_conditions = conditions(ismember(conditions(:,1), [1,2,3]) & conditions(:,2)==-90,:);
                    catch
                        keyboard
                    end
                end
                [~, condType] = ismember(conditions, use_conditions, 'rows');

            end
            
            
            %condIdx = [17, 19, 25];
            for iCond = 1:size(use_conditions,1)
                baseline_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    3, condType == iCond, 1:100), 3), 4)) .* 100;
                resp_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    3, condType == iCond, 260:360), 3), 4)) .* 100;
                std_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    3, condType == iCond, 1:100), 3), 4)) .* 100;
                psth_all(iUnit,iCond,:) = squeeze(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
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
            
            firing_rate_per_condition{iRegion, iDataset}(iUnit,:) = resp_per_cond;
            bl_firing_rate_per_condition{iRegion, iDataset}(iUnit,:) = baseline_per_cond;
            std_per_condition{iRegion, iDataset}(iUnit,:) = std_per_cond;
        end

        subplot(3, 3, (3 * (iRegion - 1))+1)
        %subplot(131)
        smooth_filt = [1, 10]; % (units x frames)
        img1 = (squeeze(psth_all(:,1,:)) - nanmean(psth_all(:,1,1:100),3))  ./ (nanstd(psth_all(:,1,1:100),[],3) );
        [~, cell_idx] = sort(nanmean(img1(:, 260:360), 2));
        img1_sm = conv2(img1(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img1(cell_idx, :)), ones(smooth_filt), 'same');
        imagesc(img1_sm)
        caxis([-4 4])

        subplot(3, 3, (3 * (iRegion - 1))+3)
        img2 = (squeeze(psth_all(:,2,:)) - nanmean(psth_all(:,2,1:100),3))  ./ (nanstd(psth_all(:,2,1:100),[],3) );
        img2_sm = conv2(img2(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img2(cell_idx, :)), ones(smooth_filt), 'same');
        imagesc(img2_sm)
        caxis([-4 4])

        subplot(3, 3, (3 * (iRegion - 1))+2)
        img3 = (squeeze(psth_all(:,3,:)) - nanmean(psth_all(:,3,1:100),3))  ./ (nanstd(psth_all(:,3,1:100),[],3) );
        img3_sm = conv2(img3(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(img3(cell_idx, :)), ones(smooth_filt), 'same');
        imagesc(img3_sm)
        caxis([-4 4])

        colormap(brewermap([],'*RdBu'))
        

        psth_per_cond{iRegion, iDataset}(:,:,1) = img1_sm;
        psth_per_cond{iRegion, iDataset}(:,:,2) = img2_sm;
        psth_per_cond{iRegion, iDataset}(:,:,3) = img3_sm;
        
    end
end

%% correlation 
dot_size = 3;
for iDataset = 3:4

end

%% firing rate
dot_size = 3;

for iDataset = 3:4
    figure();
    for iRegion = 1:3
        % remove zero rows 
        nonEmpty_rows = sum(firing_rate_per_condition{iRegion, iDataset}, 2)>0;
        % plot
        subplot(3,3,(iRegion-1)*3+1)
        act_go2 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,3)- bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,3))./...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows,3);%go2
        act_stim3 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,2)- bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,2))./...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows,2);%nogo
        act_go1 = (firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,1)- bl_firing_rate_per_condition{iRegion, iDataset}(nonEmpty_rows,1))./...
            std_per_condition{iRegion, iDataset}(nonEmpty_rows,1);%go1


        scatter(act_go2, act_go1, dot_size , 'filled')
        xlim([0.001, 100])
        ylim([0.001, 100])
        line([0.001, 100],[0.001, 100])
        set(gca, 'XScale', 'log', 'YScale', 'log');

        subplot(3,3,(iRegion-1)*3+2)
        scatter(act_go2, act_stim3, dot_size , 'filled')
        xlim([0.001, 100])
        ylim([0.001, 100])
        line([0.001, 100],[0.001, 100])
        set(gca, 'XScale', 'log', 'YScale', 'log');

        subplot(3,3,(iRegion-1)*3+3)
        scatter(act_go1, act_stim3, dot_size , 'filled')
        xlim([0.001, 100])
        ylim([0.001, 100])
        line([0.001, 100],[0.001, 100])
        set(gca, 'XScale', 'log', 'YScale', 'log');
    end
    
  %  prettify_plot('YLimits','keep','XLimits','keep')
end


