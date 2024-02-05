
%% cl_selectivity

% from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8114876/:
% Might a calculation that is more robust to trial-to-trial variability reduce
% the sensitivity of measurements to inclusion criteria or CV? We recalculated
% OSI, DSI and TF with cross-validation, using half of the trials to identify
% the stimulus condition that evoked the largest mean responses (grating
% direction and TF) and then calculated OSI, DSI and TF for these preferred
% conditions from the other half of the trials. The overall effect of
% including more neurons based on their CV on the cross-validated metrics
% across different areas was similar to that on the non-cross-validated
% metrics (Fig. 5). The notable difference is that the noisy neurons in
% the lowest decile of robustness no longer have high DSI or OSI values,
% but are shifted to much lower values (Fig. 5F,J). This difference is
% also reflected in the fact that the overall curves are shifted to lower
% values (compare Figs. 5E,I and Fig. 4E,I). Thus, while more statistically
% robust metrics calculated through cross-validation likely better reflect
% the true values of the population, they do not reduce the impact of
% selection on those metrics.

% from hikosaka https://www.jneurosci.org/content/39/9/1709:
% Object selectivity of visual neurons in the caudal putamen and caudate tail.
% A, In the passive-viewing task, fractal objects (no value associated) were
% presented sequentially in the preferred location of the neuron, while the
% monkey was fixating at the center. The reward was not associated with
% particular objects. B, Object-selective responses of a PUTt neuron from
% monkey S to fractal objects during the passive-viewing task. C–E, Distributions
% of the object selectivity index among neurons in cdPUT (C), PUTt (D), and CDt
% (E). The neurons showing significant object selectivity are indicated by black
% (ANOVA, p < 0.05). F, Object selectivity in each region (mean and 95%
% confidence interval). The across-region difference is also shown (Scheffé
% test, *p < 0.05).

% for each unit, get maximum response on half trials, and then calculate
% selectivity index based on this
% matrix of neurons mean response to each condition. use absolute value
% (mean subtracted)

% QQ homogenize so same number conds 


datasetlocations = {'/home/julie/Dropbox/MATLAB/naive_data1.mat', ...%gratings
    '/home/julie/Dropbox/MATLAB/naive_data1.mat',...
    '/home/julie/Dropbox/MATLAB/naive_data2.mat', ...%ori
    '/home/julie/Dropbox/MATLAB/naive_data3.mat', ...%nat images
    '/home/julie/Dropbox/MATLAB/naive_data6.mat', ...%nat images
    '/home/julie/Dropbox/MATLAB/naive_data4.mat', ...%cw stims
    '/home/julie/Dropbox/MATLAB/naive_data5.mat', ...%cw stims
    '/home/julie/Dropbox/MATLAB/gogogo_data2.mat', ...%cw stims go go go
    '/home/julie/Dropbox/MATLAB/gonogo_data2.mat'};%cw stims go no go
conditionsIndex = [1, 2, 6, 3, 5, 4];
selectivity_index = cell(6,3);
for iDataset = [1:3]

    % for each condition, measure
    passive_data_per_cond = load(datasetlocations{iDataset});
    keep passive_data_per_cond regions conditionsIndex iDataset datasetlocations selectivity_index
    n_conditions = 1;
    %conditionsIndex = 6;%[1,2,6];
    for iRegion = [1, 2, 3]

        %for iCondition = 1:n_conditions %QQ for 4th, combine 5th, + only use central ?
        [regionClassification, unitClassification] = cl_subsection_region(passive_data_per_cond.unit_coords, ...
            passive_data_per_cond.unit_area, passive_data_per_cond.pss, passive_data_per_cond.templateDuration, passive_data_per_cond.propLongISI);
        unique_subReg = unique(regionClassification(passive_data_per_cond.unit_area == iRegion));
        unique_subReg(unique_subReg == "") = [];

        if iRegion == 1
            curr_units = find(passive_data_per_cond.unit_area == iRegion & (passive_data_per_cond.unitType' == 1 | passive_data_per_cond.unitType' == 2) & ...
                passive_data_per_cond.pvalue_shuffled_005{1, conditionsIndex(iDataset)}' == 1 & unitClassification == 'MSN' & regionClassification == 'dorsomedial_striatum');
        else
            curr_units = find(passive_data_per_cond.unit_area == iRegion & (passive_data_per_cond.unitType' == 1 | passive_data_per_cond.unitType' == 2) & ...
                passive_data_per_cond.pvalue_shuffled_005{1, conditionsIndex(iDataset)}' == 1);
        end

       % selectivity_index{iDataset,iRegion} = nan(size(curr_units, 1), 1);

        use_conditions = 1:size(passive_data_per_cond.psth{conditionsIndex(iDataset)}, 3);
        nonEmptyRecs = find(~cellfun(@isempty, passive_data_per_cond.psth_conditions_all));
        unique_recs = unique(passive_data_per_cond.animal_thisDate_site_shank, 'rows');

        for iUnit = 1:size(curr_units, 1)

            thisUnit = curr_units(iUnit);
            [~, iRec]= ismember(passive_data_per_cond.animal_thisDate_site_shank(thisUnit,:), unique_recs, 'rows');

            conditions = passive_data_per_cond.psth_conditions_all{nonEmptyRecs(iRec)};

            if iDataset == 1 % gratings
                conditions = conditions(:,2);
                use_conditions = unique(conditions);
                [~, condType] = ismember(conditions, use_conditions);
            elseif iDataset == 2 % orientations
                 conditions = conditions(:,3);
                use_conditions = unique(conditions);
                [~, condType] = ismember(conditions, use_conditions);
            elseif iDataset == 3 % locations
                conditions = conditions(:,1);
                use_conditions = unique(conditions);
                [~, condType] = ismember(conditions, use_conditions);
            elseif iDataset == 4 % nat images 
                conditions = conditions(:,1);
                use_conditions = unique(conditions);
                [~, condType] = ismember(conditions, use_conditions);
            elseif iDataset == 5 % nat images
                conditions = conditions(:,1);
                use_conditions = unique(conditions);
                [~, condType] = ismember(conditions, use_conditions);
            elseif iDataset == 6 % cw stims
            elseif iDataset == 7 % cw stims
            elseif iDataset == 8 || iDataset == 9 % cw stims in tasks 
            end

            
            for iCond = 1:size(use_conditions,1)
                baseline_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    1, condType == iCond, 1:200), 3), 4)) .* 100;
                baseline_sub_average_per_cond(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    1, condType == iCond, 250:450), 3), 4)) .* 100;
    
                max_half_trials = find(baseline_sub_average_per_cond == max(baseline_sub_average_per_cond));
                if length(max_half_trials) > 1
                    max_half_trials = max_half_trials(1);
                end
    
                baseline_per_cond_cv(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    2, condType == iCond, 1:200), 3), 4)) .* 100;
                baseline_sub_average_per_cond_cv(iCond) = squeeze(nanmean(nanmean(passive_data_per_cond.psth{conditionsIndex(iDataset)}(thisUnit,...
                    2, condType == iCond, 250:450), 3), 4)) .* 100;
            end
            
            selectivity_index{iDataset, iRegion}(iUnit) = abs((baseline_sub_average_per_cond_cv(max_half_trials) - ...
                nanmean(baseline_sub_average_per_cond_cv))./ ...
                nanmax(baseline_sub_average_per_cond_cv)); % (c.v. max  - mean ) / max



        end

    end
end


cl_plottingSettings;

figure();
datatype_sets = [1, 1; 2, 2; 3, 3; 4, 5; 6, 7; 8, 8; 9, 9];
for iDatatype = 1:3
    for iRegion = 1:size(regions, 2)
        n_conditions = length(unique(datatype_sets(iDatatype,:)));
        datasets = datatype_sets(iDatatype,:);
        
        if n_conditions ==1 
            this_selec_idx = selectivity_index{datasets(1), iRegion}(:);
        elseif n_conditions ==2
            this_selec_idx = [selectivity_index{datasets(1), iRegion}(:); selectivity_index{datasets(2), iRegion}(:)];
        end

            subplot( 1, 4,iDatatype)
            hold on;
            [N, edges, bin] = histcounts(this_selec_idx, 0:0.05:1);
            unit_n = sum(~isnan(this_selec_idx));
            hist_bin = edges(1:end-1) + diff(edges) ./ 2;
            stairs(hist_bin, N./unit_n, 'Color', regionColors{iRegion}, 'LineWidth', 2, 'LineStyle', regionLineStyle{iRegion})

            if iRegion == 1 && iDatatype == 1
                ylabel('fraction of cells')
                xlabel('c.v. selectivity index')
            end
            axis square;
            % add median
            ylim([0, 0.4])
            sem = nanstd(this_selec_idx) ./ sqrt(length(this_selec_idx(~isnan(this_selec_idx))));
            line([nanmedian(this_selec_idx) - sem, ...
                nanmedian(this_selec_idx) + sem], ...
                [0.25 + 0.05 * iRegion, 0.25 + 0.05 * iRegion], 'Color', regionColors{iRegion}, 'LineWidth', 2, 'LineStyle', regionLineStyle{iRegion});
            scatter(nanmedian(this_selec_idx), ...
                0.25+0.05*iRegion, 25, regionColors{iRegion});

            makepretty;
            

            index_concat{iRegion}=this_selec_idx;
       

    end

    %  data12 = table;
    % data12.selectivityIndex = [selectivity_index{1}(:); selectivity_index{2}(:)];
    % data12.brainRegion = [ones(length(selectivity_index{1}), 1); ones(length(selectivity_index{2}), 1) .* 2];
    % lme12 = fitlme(data12, 'selectivityIndex ~ brainRegion');
    % 
    % data13 = table;
    % data13.selectivityIndex = [selectivity_index{1}(:); selectivity_index{3}(:)];
    % data13.brainRegion = [ones(length(selectivity_index{1}), 1); ones(length(selectivity_index{3}), 1) .* 2];
    % lme13 = fitlme(data13, 'selectivityIndex ~ brainRegion');
    % 
    % data23 = table;
    % data23.selectivityIndex = [selectivity_index{2}(:); selectivity_index{3}(:)];
    % data23.brainRegion = [ones(length(selectivity_index{2}), 1) .* 1; ones(length(selectivity_index{3}), 1) .* 2];
    % lme23 = fitlme(data23, 'selectivityIndex ~ brainRegion');

    data = table;
    data.selectivityIndex = [index_concat{1}(:); index_concat{2}(:); index_concat{3}(:)];
    data.brainRegion = [ones(length(index_concat{1}), 1); ones(length(index_concat{2}), 1) .* 2; ones(length(index_concat{3}), 1) .* 3];
    lme = fitlme(data, 'selectivityIndex ~ brainRegion');
    disp(lme)
    %emm = emmeans(lme, {'brainRegion'});
    [p, t, stats] = anova1(data.selectivityIndex, data.brainRegion);
    [c, m, h, gnames] = multcompare(stats);
end

