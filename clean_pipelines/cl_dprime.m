function [d_prime, d_prime_session_num, d_prime_animal_num, d_prime_session_fraction] = cl_dprime(task_data, idx, keepVis, keepUnits, plot_regions, plotMe)

%plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data.av_per_trial, 2)
    session_nCells(iSession) = size(task_data.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);

clearvars d_prime pooled_sd_all ci
thesePairs = [1, 2; 1, 3; 2, 3];
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
            if size(for_baseline_per_neuron_per_condition, 2) == 39
                cond_inds = [10, 34; 16, 34; 16, 10];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 26
                cond_inds = [7, 23; 11, 23; 11, 7];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 4
                cond_inds = [2, 4; 3, 4; 3, 2];
            else
                disp('wtf')
            end
            clearvars d_prime_session cond_fr keepIdx
            
            if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron > 0)))
               
                
                if ~isempty(activity_per_trial_neuron)
                if size(activity_per_trial_neuron, 2) > 2 && size(activity_per_trial_neuron, 1) > 2 
                    for iNeuron = 1:size(activity_per_trial_neuron, 2)
                        for iPair = 1:3
                            trials_1 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair, 1), :), 'rows');
                            trials_2 = ismember(theseTrialTypes, trialTypes(thesePairs(iPair, 2), :), 'rows');

                            baseline_1 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair, 1), :)./0.001);

                            baseline_2 = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair, 2), :)./0.001);
                            if baseline_1 > 0.1
                                average_stim_1 = (nanmean(activity_per_trial_neuron(trials_1, iNeuron)./0.001)); % QQ baseline - normalize ?
                                average_stim_2 = (nanmean(activity_per_trial_neuron(trials_2, iNeuron)./0.001));

                                sd_stim_1 = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001); %./sqrt(sum(trials_1));
                                sd_stim_2 = nanstd(activity_per_trial_neuron(trials_2, iNeuron)./0.001); %./sqrt(sum(trials_2));
                                n_stim_1 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                                n_stim_2 = length(activity_per_trial_neuron(trials_1, iNeuron))-1;
                                %pooled_sd = sqrt((sd_stim_1 * sd_stim_1 + sd_stim_2 * sd_stim_2)./2) + 1;
                                pooled_sd = sqrt((n_stim_1*sd_stim_1 * sd_stim_1 + n_stim_2*sd_stim_2 * sd_stim_2)./(n_stim_1+n_stim_2)) + 0.01;
                                pooled_sd_all{iRegion}(iNeuron + unitCount, iPair) = pooled_sd;
                                d_prime{iRegion}(iNeuron + unitCount, iPair) = abs(average_stim_1-average_stim_2) / pooled_sd;
                                d_prime_session{iRegion}(iNeuron, iPair) = abs(average_stim_1-average_stim_2) / pooled_sd;
                                ci{iRegion}(iNeuron + unitCount, iPair) = (average_stim_1 - average_stim_2) / (average_stim_1 + average_stim_2 + 0.1);
                                d_prime_session_num{iRegion}(iNeuron + unitCount, iPair) = iSession;
                                tempDur{iRegion}(iNeuron + unitCount, iPair) = task_data.wvDur(these_units_session_all(iNeuron));
                                pss{iRegion}(iNeuron + unitCount, iPair) = task_data.pss(these_units_session_all(iNeuron));
                                propISI{iRegion}(iNeuron + unitCount, iPair) = task_data.propISI(these_units_session_all(iNeuron));
                                if d_prime{iRegion}(iNeuron + unitCount, iPair) > 1.3
                                    %keyboard;
                                    % figure();
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(1,1), :)), 'gaussian', [0, 50])); hold on;
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(1,2), :)), 'gaussian', [0, 50]));
                                    % plot(smoothdata(squeeze(av_psth_here(iNeuron, cond_inds(2,1), :)), 'gaussian', [0, 50]));
                                    % title([num2str(abs(average_stim_1-average_stim_2)),', ', num2str(pooled_sd)])
                                end
                                %qq thissession doesn't work. wierd! 
                                %d_prime_z(iNeuron + unitCount, iPair, iRegion) = zscore((activity_per_trial_neuron(trials_1, iNeuron)))...
                                %    - ztrans((activity_per_trial_neuron(trials_2, iNeuron)));

                            else
                                d_prime_session{iRegion}(iNeuron, iPair) = NaN;
                                d_prime{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                d_prime_session_num{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                tempdur{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                pss{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                propISI{iRegion}(iNeuron + unitCount, iPair) = NaN;
                               
                            end
                        end
                        
                    end
                    unitCount = unitCount + size(activity_per_trial_neuron, 2);
                else
                end
                end
                for iPair = 1:3
                    try
                        d_prime_session_fraction(iRegion, iSession, iPair) = sum(d_prime_session{iRegion}(:, iPair) > 0.5) ./ size(d_prime_session{iRegion}, 1);
                        d_prime_session_fraction_median(iRegion, iSession, iPair) = nanmean(d_prime_session{iRegion}(:, iPair));
                    catch
                        d_prime_session_fraction(iRegion, iSession, iPair) = NaN;
                        d_prime_session_fraction_median(iRegion, iSession, iPair) = NaN;

                    end
                end
                d_prime_animal_num{iRegion, iSession} = unique_combs(iSession,1);

                % figure();
                % for iPair=1:3
                %    subplot(3,1,iPair)
                %    histogram(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair),50) % 1.4142   -1.4142 = one of stims has 0 spikes
                %    nanmedian(d_prime(unitCount+1:unitCount + size(activity_per_trial_neuron,2),iPair))
                %    title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2))])
                % end

            else
                d_prime_session_fraction(iRegion, iSession, 1:3) = NaN;
                d_prime_session_fraction_median(iRegion, iSession, 1:3) = NaN;
            end

        else
            d_prime_session_fraction(iRegion, iSession, 1:3) = NaN;
            d_prime_session_fraction_median(iRegion, iSession, 1:3) = NaN;
        end

    end

end

if plotMe
% plot d prime
figure();
for iRegion = 1:size(plot_regions, 2)
    for iPair = 1:3
        subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)))
        % if 0, remove
        %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
        kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))));

        histogram(d_prime{iRegion}(kp, iPair), [0:0.05:2]) % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(d_prime(kp,iPair,iRegion))
        %nanmean(d_prime(kp,iPair,iRegion))
        title([num2str(thesePairs(iPair, 1)), ' vs ', num2str(thesePairs(iPair, 2)), ', mean= ', num2str(nanmean(abs(d_prime{iRegion}(kp, iPair))))])
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
        for iPair = 1:3
            subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)))
            % if 0, remove
            %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
            kp = find(abs(ci{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(ci{iRegion}(:, iPair))));

            histogram(ci{iRegion}(kp, iPair), 50) % 1.4142   -1.4142 = one of stims has 0 spikes
            %nanmedian(d_prime(kp,iPair,iRegion))
            %nanmean(d_prime(kp,iPair,iRegion))
            title([num2str(thesePairs(iPair, 1)), ' vs ', num2str(thesePairs(iPair, 2)), ', mean= ', num2str(nanmean(abs(ci{iRegion}(kp, iPair))))])
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
        % kp = find(abs(d_prime(:,iPair,iRegion)) < 1.4141 & abs(d_prime(:,iPair,iRegion)) ~= 0);
        %  kp = find( abs(d_prime(:,iPair,iRegion)) ~= 0);

        scatter(squeeze(d_prime{iRegion}(:, 1)), squeeze(d_prime{iRegion}(:, 2)), 4, 'filled');
        hold on; % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(d_prime(kp,iPair,iRegion))
        %nanmean(d_prime(kp,iPair,iRegion))
        %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime(kp,iPair,iRegion))))])
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
        % kp = find(abs(d_prime(:,iPair,iRegion)) < 1.4141 & abs(d_prime(:,iPair,iRegion)) ~= 0);
        %  kp = find( abs(d_prime(:,iPair,iRegion)) ~= 0);


        scatterHistDiff(squeeze(d_prime{iRegion}(:, 1)), squeeze(d_prime{iRegion}(:, 2)), '', '', rgb('Blue'), 0.5)
        %scatter(squeeze(d_prime{iRegion}(:,1)),squeeze(d_prime{iRegion}(:,2)),4,'filled'); hold on; % 1.4142   -1.4142 = one of stims has 0 spikes
        %nanmedian(d_prime(kp,iPair,iRegion))
        %nanmean(d_prime(kp,iPair,iRegion))
        %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime(kp,iPair,iRegion))))])
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