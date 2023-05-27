
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

% for each unit, get maximum response on half trials, and then calculate
% selectivity index based on this
% matrix of neurons mean response to each condition. use absolute value
% (mean subtracted)

% for each condition, measure
keep passive_data regions
selectivity_index = cell(size(regions,2),1);
n_conditions = 4;
for iRegion = 1:size(regions,2)
    curr_units = find(passive_data.unit_area == iRegion);
    selectivity_index{iRegion} = nan(size(curr_units, 1), 4);

    for iCondition = 1:n_conditions %QQ for 4th, combine 5th, + only use central ? 
        if iCondition == 4 
            use_conditions = 2:2:26;
        else
            use_conditions = 1:size(passive_data.psth{iCondition},3);
        end
        for iUnit = 1:size(curr_units, 1)
            thisUnit = curr_units(iUnit);
            baseline_per_cond = squeeze(nanmean(passive_data.psth{iCondition}(thisUnit, 1, use_conditions, 1:50), 4));
            baseline_sub_average_per_cond = arrayfun(@(x) nanmean(abs(squeeze(passive_data.psth{iCondition}(thisUnit, 1, x, 55:70)) ...
                -baseline_per_cond(x))), 1:size(baseline_per_cond, 1));

            max_half_trials = find(baseline_sub_average_per_cond == max(baseline_sub_average_per_cond));
            if length(max_half_trials) > 1
                max_half_trials = max_half_trials(1);
            end

            baseline_per_cond_cv = squeeze(nanmean(passive_data.psth{iCondition}(thisUnit, 2, use_conditions, 1:50), 4));
            baseline_sub_average_per_cond_cv = arrayfun(@(x) nanmean(abs(squeeze(passive_data.psth{iCondition}(thisUnit, 2, x, 55:70)) ...
                -baseline_per_cond_cv(x))+baseline_per_cond_cv(x)), 1:size(baseline_per_cond_cv, 1));
            selectivity_index{iRegion}(thisUnit, iCondition) = abs((baseline_sub_average_per_cond_cv(max_half_trials) - nanmean(baseline_sub_average_per_cond_cv))./ ...
                nanmax(baseline_sub_average_per_cond_cv)); % (c.v. max  - mean ) / max

        end
    end
end

theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

figure();
for iRegion = 1:size(regions,2)
    for iCondition = 1:n_conditions
        subplot(1, n_conditions, iCondition)
        hold on;
        [N,edges,bin] = histcounts(selectivity_index{iRegion}(:,iCondition), 0:0.05:1);
        hist_bin = edges(1:end-1) + diff(edges) ./ 2;
        stairs(hist_bin, N./length(bin), 'Color', theseColors{iRegion}, 'LineWidth', 2)
        
        if iRegion == 1 && iCondition == 1
            ylabel('fraction of cells')
            xlabel('c.v. selectivity index')
        end
        axis square;
        makepretty; 
    end
end
