function [vis_resp, vis_resp_session_num, vis_resp_animal_num, vis_resp_session_fraction, vis_resp_full] = cl_visResp(task_data, idx, keepVis, keepUnits, plot_regions, plotMe)

%plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data.av_per_trial, 2)
    session_nCells(iSession) = size(task_data.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);


theseStims = [1; 2; 3];
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
            if size(task_data.av_psth{idx, thisSession},3) == 800
                for_baseline_per_neuron_per_condition = task_data.av_psth{idx, thisSession}(these_units_session, :, 1:200);
            else
                for_baseline_per_neuron_per_condition = task_data.av_psth{idx, thisSession}(these_units_session, :, 1:50);
            end

            if size(for_baseline_per_neuron_per_condition, 2) == 39
                cond_inds = [10; 34; 16];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 26
                cond_inds = [7; 23; 11];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 4
                cond_inds = [2; 4; 3];
            elseif size(for_baseline_per_neuron_per_condition, 2) == 24
                disp('error')
                continue;
            elseif size(for_baseline_per_neuron_per_condition, 2) == 66
                cond_inds = [10; 34; 16];
            else
                disp('wtf')
            end
            clearvars d_prime_session cond_fr keepIdx
            
            if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron > 0)))
               
                
                if ~isempty(activity_per_trial_neuron)
                if size(activity_per_trial_neuron, 2) > 2 && size(activity_per_trial_neuron, 1) > 2 
                    for iNeuron = 1:size(activity_per_trial_neuron, 2)
                        for iPair = 1:3
                            trials_1 = ismember(theseTrialTypes, trialTypes(theseStims(iPair, 1), :), 'rows');
                           

                            baseline = nanmean(for_baseline_per_neuron_per_condition(iNeuron, ...
                                cond_inds(iPair, 1), :)./0.001);

                            
                           

                            if baseline > 0.1
                                average_stim = (nanmean(activity_per_trial_neuron(logical(trials_1), iNeuron)./0.001)); % QQ baseline - normalize ?

                                sd_stim = nanstd(activity_per_trial_neuron(trials_1, iNeuron)./0.001); %./sqrt(sum(trials_1));
                                
                                
                                %vis_resp{iRegion}(iNeuron + unitCount, iPair) = (average_stim - baseline  ) / abs(average_stim + baseline+ 0.1);
                                %vis_resp_session{iRegion}(iNeuron, iPair) = (average_stim - baseline) / abs(average_stim+baseline+ 0.1);
                                train_image = (train_image - nanmean(train_image(:, val_t_1), 2)) ./ ...
            nanstd(train_image(:, val_t_1), [], 2);
        [~, cell_idx] = sort(nanmean(train_image(:, val_t_2), 2));

                                
                                vis_resp_full{iRegion}(iNeuron + unitCount, iPair,:) = av_psth_here(iNeuron, cond_inds(iPair, 1), :);
                                %vis_resp{iRegion}(iNeuron + unitCount, iPair) = (average_stim - baseline  ) / abs(sd_stim+ 0.001);
                                %vis_resp_session{iRegion}(iNeuron, iPair) = (average_stim - baseline) / abs(sd_stim+ 0.001);
                                
                                vis_resp_session_num{iRegion}(iNeuron + unitCount, iPair) = iSession;
                                
                                ci{iRegion}(iNeuron + unitCount, iPair) = (average_stim - baseline) / (average_stim + baseline + 0.1);
                                
                                tempDur{iRegion}(iNeuron + unitCount, iPair) = task_data.wvDur(these_units_session_all(iNeuron));
                                pss{iRegion}(iNeuron + unitCount, iPair) = task_data.pss(these_units_session_all(iNeuron));
                                propISI{iRegion}(iNeuron + unitCount, iPair) = task_data.propISI(these_units_session_all(iNeuron));
                                
                                %qq thissession doesn't work. wierd! 
                                %d_prime_z(iNeuron + unitCount, iPair, iRegion) = zscore((activity_per_trial_neuron(trials_1, iNeuron)))...
                                %    - ztrans((activity_per_trial_neuron(trials_2, iNeuron)));

                            else
                                vis_resp_session{iRegion}(iNeuron, iPair) = NaN;
                                vis_resp{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                vis_resp_session_num{iRegion}(iNeuron + unitCount, iPair) = NaN;
                                vis_resp_full{iRegion}(iNeuron + unitCount, iPair,:) = nan(800,1);
                                
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
                        vis_resp_session_fraction(iRegion, iSession, iPair) = sum(vis_resp_session{iRegion}(:, iPair) > 0.5) ./ size(vis_resp_session{iRegion}, 1);
                    catch
                        vis_resp_session_fraction(iRegion, iSession, iPair) = NaN;

                    end
                end
                vis_resp_animal_num{iRegion, iSession} = unique_combs(iSession,1);

              

            else
                vis_resp_session_fraction(iRegion, iSession, 1:3) = NaN;
              
            end

        else
            vis_resp_session_fraction(iRegion, iSession, 1:3) = NaN;
            
        end

    end

end

if plotMe

end