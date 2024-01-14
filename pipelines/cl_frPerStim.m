function [fr, fr_session_num, fr_animal_num, fr_session_fraction] = cl_frPerStim(task_data, idx, keepVis, keepUnits, plot_regions, plotMe)

%plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data.av_per_trial, 2)
    session_nCells(iSession) = size(task_data.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);

clearvars fr pooled_sd_all ci
thesePairs = [1, 2, 3];
task_data.animal_day_site_shank(isnan(task_data.animal_day_site_shank)) =0;
unique_combs = unique(task_data.animal_day_site_shank, 'rows');
unitCount = 0;

for iRegion = 1:size(plot_regions, 2)
    if keepVis
        these_units = task_data.unit_area == iRegion & ...
            ismember(task_data.unitType' , keepUnits)& ...
        task_data.pvalue_shuffled_005{1, idx}' == 1;%any(task_data.pvalue{1,idx} < 0.01,2);
    else
        these_units = task_data.unit_area == iRegion & ...
            ismember(task_data.unitType' , keepUnits); % & ...
        %task_data.pvalue_shuffled_005{1,idx}' == 1;
    end
    if iRegion ==1
    end
    unitCount = 0;
    for iSession = 1:size(nonZero_idx, 2) - 1
        thisSession = nonZero_idx(iSession+1) - 1;
        trials_no_move = task_data.trial_types{idx, thisSession}(task_data.no_move_trials{idx, thisSession}, :);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));


        if ~isempty(trials_no_move)
            theseTrials = ismember(trials_no_move, trialTypes, 'rows'); %QQ
            theseTrialTypes = trials_no_move(theseTrials, :);
            these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
            these_units_session_all = find(these_units);
            these_units_session_all = these_units_session_all(...
                these_units_session_all >= session_cumCells(thisSession) & these_units_session_all <= session_cumCells(thisSession+1)-1);
            %figure();
            %i=2;
            %imagesc(squeeze(task_data.av_psth{idx,thisSession}(:,i,:)))
            activity_per_trial_neuron = task_data.av_per_trial{idx, thisSession}(these_units_session, theseTrials)';
            av_psth_here = task_data.av_psth{idx, thisSession}(these_units_session, :, :);
            for_baseline_per_neuron_per_condition = task_data.av_psth{idx, thisSession}(these_units_session, :, 1:200);
            %disp(size(for_baseline_per_neuron_per_condition))
            if size(for_baseline_per_neuron_per_condition, 2) == 39 || size(for_baseline_per_neuron_per_condition, 2) == 66
                cond_inds = [10, 34, 16];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 26
                cond_inds = [7, 23, 11];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 4
                cond_inds = [2, 4, 3];
            else
                disp('wtf')
                continue;
            end
            clearvars fr_session cond_fr keepIdx
            
            if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron > 0)))
               
                
                if ~isempty(activity_per_trial_neuron)
                if size(activity_per_trial_neuron, 2) > 2 && size(activity_per_trial_neuron, 1) > 2 
                    for iNeuron = 1:size(activity_per_trial_neuron, 2)
                        for iStim = 1:3
                            trials_1 = ismember(theseTrialTypes, trialTypes(thesePairs(iStim), :), 'rows');
                            %trials_2 = ismember(theseTrialTypes, trialTypes(thesePairs(iStim, 2), :), 'rows');

                            baseline_1 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iStim), :)./0.001);

                            %baseline_2 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                            %    cond_inds(iStim, 2), :)./0.001);
                            if baseline_1 > 0.1
                                average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ?
                               % average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001));

                                sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001); %./sqrt(sum(trials_1));
                                %sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001); %./sqrt(sum(trials_2));
                                n_stim_1 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                               % n_stim_2 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                                %pooled_sd = sqrt((sd_stim_1 * sd_stim_1 + sd_stim_2 * sd_stim_2)./2) + 1;
                               % pooled_sd = sqrt((n_stim_1*sd_stim_1 * sd_stim_1 + n_stim_2*sd_stim_2 * sd_stim_2)./(n_stim_1+n_stim_2)) + 0.01;
                               % pooled_sd_all{iRegion}(iNeuron + unitCount, iStim) = pooled_sd;
                                
                                fr{iRegion}(iNeuron + unitCount, iStim) = (average_stim_1 - baseline_1)./sd_stim_1;
                                fr_session{iRegion}(iNeuron, iStim) = (average_stim_1 - baseline_1)./sd_stim_1;
                                
                                
                                %ci{iRegion}(iNeuron + unitCount, iStim) = (average_stim_1 - average_stim_2) / (average_stim_1 + average_stim_2 + 0.1);
                                fr_session_num{iRegion}(iNeuron + unitCount, iStim) = iSession;
                                try
                                tempDur{iRegion}(iNeuron + unitCount, iStim) = task_data.templateDuration(these_units_session_all(iNeuron));
                                catch
                                    tempDur{iRegion}(iNeuron + unitCount, iStim) = task_data.wvDur(these_units_session_all(iNeuron));
                                end
                                pss{iRegion}(iNeuron + unitCount, iStim) = task_data.pss(these_units_session_all(iNeuron));
                                try
                                propISI{iRegion}(iNeuron + unitCount, iStim) = task_data.propLongISI(these_units_session_all(iNeuron));
                                catch
                                    propISI{iRegion}(iNeuron + unitCount, iStim) = task_data.propISI(these_units_session_all(iNeuron));
                                end
                                if fr{iRegion}(iNeuron + unitCount, iStim) > 1.3
                                    %keyboard;
                                    % figure();
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(1,1), :)), 'gaussian', [0, 50])); hold on;
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(1,2), :)), 'gaussian', [0, 50]));
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(2,1), :)), 'gaussian', [0, 50]));
                                    % title([num2str(abs(average_stim_1-average_stim_2)),', ', num2str(pooled_sd)])
                                end
                                %qq thissession doesn't work. wierd! 
                                %fr_z(iNeuron + unitCount, iPair, iRegion) = zscore((activity_per_trial_neuron(trials_1, iNeuron)))...
                                %    - ztrans((activity_per_trial_neuron(trials_2, iNeuron)));

                            else
                                fr_session{iRegion}(iNeuron, iStim) = NaN;
                                fr{iRegion}(iNeuron + unitCount, iStim) = NaN;
                                fr_session_num{iRegion}(iNeuron + unitCount, iStim) = NaN;
                                tempdur{iRegion}(iNeuron + unitCount, iStim) = NaN;
                                pss{iRegion}(iNeuron + unitCount, iStim) = NaN;
                                propISI{iRegion}(iNeuron + unitCount, iStim) = NaN;
                               
                            end
                        end
                        
                    end
                    unitCount = unitCount + size(activity_per_trial_neuron, 2);
                else
                end
                end
                for iStim = 1:3
                    try
                        fr_session_fraction(iRegion, iSession, iStim) = sum(fr_session{iRegion}(:, iStim) > 0.5) ./ size(fr_session{iRegion}, 1);
                        fr_session_fraction_median(iRegion, iSession, iStim) = nanmean(fr_session{iRegion}(:, iStim));
                    catch
                        fr_session_fraction(iRegion, iSession, iStim) = NaN;
                        fr_session_fraction_median(iRegion, iSession, iStim) = NaN;

                    end
                end
                fr_animal_num{iRegion, iSession} = unique_combs(iSession,1);

                % figure();
                % for iPair=1:3
                %    subplot(3,1,iPair)
                %    histogram(fr(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes
                %    nanmedian(fr(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair))
                %    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2))])
                % end

            else
                fr_session_fraction(iRegion, iSession, 1:3) = NaN;
                fr_session_fraction_median(iRegion, iSession, 1:3) = NaN;
            end

        else
            fr_session_fraction(iRegion, iSession, 1:3) = NaN;
            fr_session_fraction_median(iRegion, iSession, 1:3) = NaN;
        end

    end

end

if plotMe
% plot d prime
figure();
for iRegion = 1:size(plot_regions, 2)
    for iStim = 1:3
        subplot(3, size(plot_regions, 2), iStim+(iRegion - 1)*(size(plot_regions, 2)))
        % if 0, remove
        %kp = find(abs(fr{iRegion}(:,iPair)) < 1.4141 & abs(fr{iRegion}(:,iPair)) ~= 0);
        kp = find(abs(fr{iRegion}(:, iStim)) ~= 0 & ~isinf(abs(fr{iRegion}(:, iStim))));

        histogram(fr{iRegion}(kp, iStim), [0:0.05:2]) % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(fr(kp,iPair,iRegion))
        %nanmean(fr(kp,iPair,iRegion))
        title([num2str(thesePairs(iStim)), ' vs ', num2str(thesePairs(iStim, 2)), ', mean= ', num2str(nanmean(abs(fr{iRegion}(kp, iStim))))])
        xlabel('d-prime')
        ylabel('# of neurons')
        xlim([0, 2])
        %makepretty_lite;
    end
end

more = 0;
if more
    % plot sel. idx
    figure();
    for iRegion = 1:size(plot_regions, 2)
        for iStim = 1:3
            subplot(3, size(plot_regions, 2), iStim+(iRegion - 1)*(size(plot_regions, 2)))
            % if 0, remove
            %kp = find(abs(fr{iRegion}(:,iPair)) < 1.4141 & abs(fr{iRegion}(:,iPair)) ~= 0);
            kp = find(abs(ci{iRegion}(:, iStim)) ~= 0 & ~isinf(abs(ci{iRegion}(:, iStim))));

            histogram(ci{iRegion}(kp, iStim), 50) % 1.4142   -1.4142 = one of stims has 0 spikes
            %nanmedian(fr(kp,iPair,iRegion))
            %nanmean(fr(kp,iPair,iRegion))
            title([num2str(thesePairs(iStim)), ' vs ', num2str(thesePairs(iStim, 2)), ', mean= ', num2str(nanmean(abs(ci{iRegion}(kp, iStim))))])
            xlabel('Sel. idx')
            ylabel('# of neurons')
            %xlim([-1.45, 1.45])
            %makepretty_lite;
        end
    end


    figure();
    for iRegion = 1:size(plot_regions, 2)
        %for iPair=1:3
        subplot(size(plot_regions, 2), 1, iRegion)
        % if 0, remove
        % kp = find(abs(fr(:,iPair,iRegion)) < 1.4141 & abs(fr(:,iPair,iRegion)) ~= 0);
        %  kp = find( abs(fr(:,iPair,iRegion)) ~= 0);

        scatter(squeeze(fr{iRegion}(:, 1)), squeeze(fr{iRegion}(:, 2)), 4, 'filled');
        hold on; % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(fr(kp,iPair,iRegion))
        %nanmean(fr(kp,iPair,iRegion))
        %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(fr(kp,iPair,iRegion))))])
        xlabel('1 vs 2')
        ylabel('1 vs 3')
        %xlim([-1.4, 1.4])
        % ylim([-1.4, 1.4])


        set(gca, 'DataAspectRatio', [1, 1, 1])
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        xlimits = xlim;
        ylimits = ylim;
        maxLim = max([xlimits, ylimits]);
        minLim = min([xlimits, ylimits]);
        line([minLim, maxLim], [minLim, maxLim], 'Color', rgb('Red'));
        %makepretty;
        %get()


        %end
    end

    figure();
    for iRegion = 1:size(plot_regions, 2)
        %for iPair=1:3
        subplot(size(plot_regions, 2), 1, iRegion)
        % if 0, remove
        % kp = find(abs(fr(:,iPair,iRegion)) < 1.4141 & abs(fr(:,iPair,iRegion)) ~= 0);
        %  kp = find( abs(fr(:,iPair,iRegion)) ~= 0);


        scatterHistDiff(squeeze(fr{iRegion}(:, 1)), squeeze(fr{iRegion}(:, 2)), '', '', rgb('Blue'), 0.5)
        %scatter(squeeze(fr{iRegion}(:,1)),squeeze(fr{iRegion}(:,2)),4,'filled'); hold on; % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(fr(kp,iPair,iRegion))
        %nanmean(fr(kp,iPair,iRegion))
        %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(fr(kp,iPair,iRegion))))])
        xlabel('1 vs 2')
        ylabel('1 vs 3')
        %xlim([-1.4, 1.4])
        % ylim([-1.4, 1.4])


        set(gca, 'DataAspectRatio', [1, 1, 1])
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        xlimits = xlim;
        ylimits = ylim;
        maxLim = max([xlimits, ylimits]);
        minLim = min([xlimits, ylimits]);
        line([minLim, maxLim], [minLim, maxLim], 'Color', rgb('Red'));
        %makepretty;
        %get()


        %end
    end

    scatterHistDiff(x, y, xeb, yeb, colors, histogramBinSize)
end
end