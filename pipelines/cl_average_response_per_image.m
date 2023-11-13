%% choiceworld 
task_data_here = load('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat');

keep task_data_here
idx = 5; %2: task, 4/5 : passive
first_dataset = [0, 0; 1, -90; 2, -90; 3, -90; 4, -90; 5, -90; 6, -90; 7, -90; 8, -90; 8, 0; 9, 0; 10, 0; 11, 0; 12, 0; 13, 0];
second_dataset = [1, -90; 1, 0; 1, 90; 2, -90; 2, 0; 2, 90; 3, -90; 3, 0; 3, 90; 4, -90; 4, 0; 4, 90; 5, -90; 5, 0; 5, 90; ...
    6, -90; 6, 0; 6, 90; 7, -90; 7, 0; 7, 90; 8, -90; 8, 0; 8, 90; 9, -90; 9, 0; 9, 90; 10, -90; 10, 0; 10, 90; ...
    11, -90; 11, 0; 11, 90; 12, -90; 12, 0; 12, 90; 13, -90; 13, 0; 13, 90];
third_dataset = [0,0; 1, -90; 1, 0;  2, -90; 2, 0;  3, -90; 3, 0;  4, -90; 4, 0; 5, -90; 5, 0;  ...
    6, -90; 6, 0;  7, -90; 7, 0;  8, -90; 8, 0;9, -90; 9, 0;  10, -90; 10, 0;  ...
    11, -90; 11, 0; 12, -90; 12, 0;  13, -90; 13, 0];


plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data_here.av_per_trial, 2)
    session_nCells(iSession) = size(task_data_here.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
trialTypes = [4, -90; 12, -90; 6, -90]; %go1, go2, no go
session_cumCells = cumsum(session_nCells);

clearvars d_prime pooled_sd_all ci
thesePairs = [1, 2; 1, 3; 2, 3];

unitCount = 0;
responses_per_stim = nan(38, 3);
responses_per_stim_std = nan(38, 3);
responses_per_stim_1 = nan(38, 3);
responses_per_stim_std_1 = nan(38, 3);
responses_per_stim_2 = nan(38, 3);
responses_per_stim_std_2 = nan(38, 3);

responses_session_base = nan(size(nonZero_idx, 2)-1, 38, 3);
responses_session = nan(size(nonZero_idx, 2)-1, 38, 3);
responses_session_base_1 = nan(size(nonZero_idx, 2)-1, 38, 3);
responses_session_1 = nan(size(nonZero_idx, 2)-1, 38, 3);
responses_session_base_2 = nan(size(nonZero_idx, 2)-1, 38, 3);
responses_session_2 = nan(size(nonZero_idx, 2)-1, 38, 3);


for iRegion = 1:size(plot_regions, 2)
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' == 1 | task_data_here.unitType' == 2) & ...
         task_data_here.pvalue_shuffled_0025{1,idx}' == 1;
    unitCount = 0;
    for iSession = 1:size(nonZero_idx, 2) - 1
        thisSession = nonZero_idx(iSession+1) - 1;
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1) & ...
            nanmean(nanmean(task_data_here.av_psth{idx, thisSession},2),3)> 0.005;
        if sum(these_units_session) > 0
            % coarse way of determing contral. 
            nm_per_stim = nanmean(task_data_here.av_psth{idx, thisSession},3);
            nm_per_stim_base = nanmean(task_data_here.av_psth{idx, thisSession}(:,:,1:200),3);

            nm_per_stim_abs = abs(nanmean((nm_per_stim - nm_per_stim_base + 0.0001 ) ./ (nm_per_stim_base + 0.0001)));
            
            [~, idx_max] = max([nanmean(nm_per_stim_abs(task_data_here.psth_conditions_all{thisSession,idx}(:,2)==-90)), ...
                nanmean(nm_per_stim_abs(task_data_here.psth_conditions_all{thisSession,idx}(:,2)==90))]);
            if idx_max == 1
                contra_ids = find(task_data_here.psth_conditions_all{thisSession, idx}(:,2) == -90 | task_data_here.psth_conditions_all{thisSession, idx}(:,2) == 0);
            else
                contra_ids = find(task_data_here.psth_conditions_all{thisSession, idx}(:,2) == 90 | task_data_here.psth_conditions_all{thisSession, idx}(:,2) == 0);
          
            end
            these_units_session= find(these_units_session);
            clearvars responses_temp responses_temp_1 responses_temp_2
            for iStim = 1:size(contra_ids,1)
                
                
                stimID = task_data_here.psth_conditions_all{thisSession,idx}(contra_ids(iStim),:);
                if stimID(2) == 90
                    stimID(2) = -90;
                end
                stimRow = ismember(second_dataset, stimID, 'rows');
                for iCell = 1:size(these_units_session,1)
                    responses_temp(iCell, iStim) = abs(nanmean((task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), :) - ...
                        nanmean(task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001))));
                    responses_temp_1(iCell, iStim) = abs(nanmean((task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), :) - ...
                        nanmean(task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001))));
                    responses_temp_2(iCell, iStim) = abs(nanmean((task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), :) - ...
                        nanmean(task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), contra_ids(iStim), 1:200)+ 0.0001))));
                end
                responses_session(iSession, stimRow, iRegion) = nanmean(responses_temp(:,iStim),1);
                responses_session_1(iSession, stimRow, iRegion) = nanmean(responses_temp_1(:,iStim),1);
                responses_session_2(iSession, stimRow, iRegion) = nanmean(responses_temp_2(:,iStim),1);
                
            end
        end

    end

    responses_per_stim_std(:, iRegion) = nanstd(responses_session(:,:,iRegion), [], 1);
    responses_per_stim(:, iRegion) = nanmean(responses_session(:,:,iRegion),1);

    responses_per_stim_std_1(:, iRegion) = nanstd(responses_session_1(:,:,iRegion), [], 1);
    responses_per_stim_1(:, iRegion) = nanmean(responses_session_1(:,:,iRegion),1);

    responses_per_stim_std_2(:, iRegion) = nanstd(responses_session_2(:,:,iRegion), [], 1);
    responses_per_stim_2(:, iRegion) = nanmean(responses_session_2(:,:,iRegion),1);

 end



%% plot cross-validated PSTHs 
figure();


%% scatter intra-region
% sort based on striatum, half trials 
sort(responses_per_stim_1(:,1))
figure();
thisText = {'Str', 'GPe', 'SNr'};
for iRegion = 1:3
    subplot(1, 3, iRegion); hold on;
    scatter(responses_per_stim_1([1:3:end],iRegion), responses_per_stim_2([1:3:end],iRegion),40, 'filled')
    scatter(responses_per_stim_1([2:3:end],iRegion), responses_per_stim_2([2:3:end],iRegion), 40, rgb('MidnightBlue'), 'filled')
    %scatter(responses_per_stim_1(2:3:end,iRegion), responses_per_stim_2(2:3:end,iRegion), 'filled')
    %scatter(responses_per_stim_1(3:3:end,iRegion), responses_per_stim_2(3:3:end,iRegion), 'filled')
    
    xlim([0, 0.7])
    ylim([0, 0.7])
    
    xlabel(['increase in F.R. after stimulus onset', newline, ...
        '1/2 trials ', thisText{iRegion}])
    ylabel(['increase in F.R. after stimulus onset', newline, ...
        '1/2 trials ', thisText{iRegion}])
    %nan_infs = any(isnan(responses_per_stim_1(1:3:end, :)), 2) | any(isinf(responses_per_stim_2(1:3:end, :)), 2);
    line([0, 0.7], [0, 0.7])
    [r, p] = corrcoef(responses_per_stim_1([1:3:end, 2:3:end], 1), responses_per_stim_2([1:3:end, 2:3:end], iRegion), 'rows', 'complete');
    legend('contralateral', 'central')
prettify_plot;
end
%% scatter: inter-region increase in F.R. 
thisText = {'', 'GPe', 'SNr'};
figure();
for iRegion = 2:3
    subplot(1, 2, iRegion -1)
    hold on;
    %[~, sorted_idx] = sort(responses_per_stim(:, 1));
    scatter(responses_per_stim([1:3:end],1), responses_per_stim([1:3:end],iRegion),40 ,'filled')
    scatter(responses_per_stim([2:3:end],1), responses_per_stim([2:3:end],iRegion), 40, rgb('MidnightBlue'), 'filled')
    %scatter(responses_per_stim(2:3:end, 1), responses_per_stim(2:3:end, iRegion), 'filled') %   0*
    %scatter(responses_per_stim(3:3:end, 1), responses_per_stim(3:3:end, iRegion), 'filled') %  90*
    prettify_plot;
    xlim([0, 0.7])
    ylim([0, 0.7])
    hold on;
    xlabel(['increase in F.R. after stimulus onset', newline, ...
        'striatum'])
    ylabel(['increase in F.R. after stimulus onset', newline, ...
        thisText{iRegion}])

    nan_infs = any(isnan(responses_per_stim(:, :)), 2) | any(isinf(responses_per_stim(:, :)), 2);
    [r, p] = corrcoef(responses_per_stim([1:3:end, 2:3:end], 1), responses_per_stim([1:3:end, 2:3:end], iRegion), 'rows', 'complete');
     scatter(responses_per_stim(find(second_dataset(:,1)==4 & second_dataset(:,2)==-90), 1),...
        responses_per_stim(find(second_dataset(:,1)==4 & second_dataset(:,2)==-90), iRegion),40,rgb('IndianRed'),'filled')
     scatter(responses_per_stim(find(second_dataset(:,1)==6 & second_dataset(:,2)==-90), 1),...
         responses_per_stim(find(second_dataset(:,1)==6 & second_dataset(:,2)==-90), iRegion),40,rgb('LawnGreen'),'filled')
     scatter(responses_per_stim(find(second_dataset(:,1)==12 & second_dataset(:,2)==-90), 1),...
         responses_per_stim(find(second_dataset(:,1)==12 & second_dataset(:,2)==-90), iRegion),40,rgb('Amethyst'),'filled')

      scatter(responses_per_stim(find(second_dataset(:,1)==4 & second_dataset(:,2)==0), 1),...
        responses_per_stim(find(second_dataset(:,1)==4 & second_dataset(:,2)==0), iRegion),40,rgb('Red') ,'filled')
     scatter(responses_per_stim(find(second_dataset(:,1)==6 & second_dataset(:,2)==0), 1),...
         responses_per_stim(find(second_dataset(:,1)==6 & second_dataset(:,2)==0), iRegion),40,rgb('Green') ,'filled')
     scatter(responses_per_stim(find(second_dataset(:,1)==12 & second_dataset(:,2)==0), 1),...
         responses_per_stim(find(second_dataset(:,1)==12 & second_dataset(:,2)==0), iRegion),40,rgb('Purple') ,'filled')
   

    % errorbar(responses_per_stim(second_dataset(:,2)~=90, 1), responses_per_stim(second_dataset(:,2)==-90, iRegion),...
    %     responses_per_stim_std(second_dataset(:,2)~=90, iRegion), "LineStyle", "none", 'Color', 'k')
    % errorbar(responses_per_stim(second_dataset(:,2)~=90, 1), responses_per_stim(second_dataset(:,2)==-90, iRegion),...
    %     responses_per_stim_std(second_dataset(:,2)~=90, 1), 'horizontal', "LineStyle", "none", 'Color', 'k')
    
    line([0, 0.7], [0, 0.7])
    legend('contralateral', 'central', 'contra go1', 'contra noGo', 'contra go2', ...
        'central go1', 'central noGo', 'central go2')

end



%% nat images
indexes = [1,2,6];
textyLoads = {'/home/julie/Dropbox/MATLAB/gratings_passive_full.mat', ...
    '/home/julie/Dropbox/MATLAB/locations_passive_full.mat', '/home/julie/Dropbox/MATLAB/nat_passive_full.mat'};
for iType = 3%%:3

task_data_here = load(textyLoads{iType});

keep task_data_here iType indexes textyLoads
idx = indexes(iType); %2: task, 4/5 : passive

unique_values = unique([task_data_here.psth_conditions_all{end,idx}], 'rows');

plot_regions = [1, 2, 3]; %[1, 2, 5]; % Striatum, GPe, SNr
for iSession = 1:size(task_data_here.av_per_trial, 2)
    session_nCells(iSession) = size(task_data_here.av_per_trial{idx, iSession}, 1);
end
session_nCells = [1, session_nCells];
[~, nonZero_idx] = find(session_nCells > 0);
session_cumCells = cumsum(session_nCells);

clearvars d_prime pooled_sd_all ci

unitCount = 0;
nresponses_per_stim = nan(size(unique_values,1), 3);
nresponses_per_stim_std = nan(size(unique_values,1), 3);
nresponses_per_stim_1 = nan(size(unique_values,1), 3);
nresponses_per_stim_std_1 = nan(size(unique_values,1), 3);
nresponses_per_stim_2 = nan(size(unique_values,1), 3);
nresponses_per_stim_std_2 = nan(size(unique_values,1), 3);

responses_session_base = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);
responses_session = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);
responses_session_base_1 = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);
responses_session_1 = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);
responses_session_base_2 = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);
responses_session_2 = nan(size(nonZero_idx, 2)-1, size(unique_values,1), 3);


for iRegion = 1:size(plot_regions, 2)
    these_units = task_data_here.unit_area == iRegion & ...
        (task_data_here.unitType' == 1 | task_data_here.unitType' == 2) & ...
         task_data_here.pvalue_shuffled_0025{1,idx}' == 1;
    unitCount = 0;
    for iSession = 1:size(nonZero_idx, 2) - 1
        
        thisSession = nonZero_idx(iSession+1) - 1;
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1) & ...
            nanmean(nanmean(task_data_here.av_psth{idx, thisSession},2),3)> 0.005;
        if sum(these_units_session) > 0
            % coarse way of determing contral. 
            nm_per_stim = nanmean(task_data_here.av_psth{idx, thisSession},3);
            nm_per_stim_base = nanmean(task_data_here.av_psth{idx, thisSession}(:,:,1:200),3);

            nm_per_stim_abs = abs(nanmean((nm_per_stim - nm_per_stim_base + 0.0001 ) ./ (nm_per_stim_base + 0.0001)));
            
            %[~, idx_max] = max([nanmean(nm_per_stim_abs(task_data_here.psth_conditions_all{thisSession,idx}(:,2)==-90)), ...
            %    nanmean(nm_per_stim_abs(task_data_here.psth_conditions_all{thisSession,idx}(:,2)==90))]);
            
            these_units_session= find(these_units_session);
            clearvars responses_temp responses_temp_1 responses_temp_2
            for iStim = 1:size(task_data_here.psth_conditions_all{thisSession,idx},1)
                
                
                stimID = task_data_here.psth_conditions_all{thisSession,idx}(iStim,:);
                
                stimRow = ismember(unique_values, stimID, 'rows');
                for iCell = 1:size(these_units_session,1)
                    responses_temp(iCell, iStim) = abs(nanmean((task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), iStim, :) - ...
                        nanmean(task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001))));
                    responses_temp_1(iCell, iStim) = abs(nanmean((task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), iStim, :) - ...
                        nanmean(task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth_1{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001))));
                    responses_temp_2(iCell, iStim) = abs(nanmean((task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), iStim, :) - ...
                        nanmean(task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001)) ./...
                        (nanmean(task_data_here.av_psth_2{idx, thisSession}(these_units_session(iCell), iStim, 1:200)+ 0.0001))));
                end
                responses_session(iSession, stimRow, iRegion) = nanmean(responses_temp(:,iStim),1);
                responses_session_1(iSession, stimRow, iRegion) = nanmean(responses_temp_1(:,iStim),1);
                responses_session_2(iSession, stimRow, iRegion) = nanmean(responses_temp_2(:,iStim),1);
                
            end
        
        end

    end

    nresponses_per_stim_std(:, iRegion) = nanstd(responses_session(:,:,iRegion), [], 1);
    nresponses_per_stim(:, iRegion) = nanmean(responses_session(:,:,iRegion),1);

    nresponses_per_stim_std_1(:, iRegion) = nanstd(responses_session_1(:,:,iRegion), [], 1);
    nresponses_per_stim_1(:, iRegion) = nanmean(responses_session_1(:,:,iRegion),1);

    nresponses_per_stim_std_2(:, iRegion) = nanstd(responses_session_2(:,:,iRegion), [], 1);
    nresponses_per_stim_2(:, iRegion) = nanmean(responses_session_2(:,:,iRegion),1);

 end



%% plot cross-validated PSTHs 



%% scatter intra-region
find(nresponses_per_stim(:,1)>=0.2)

figure();
thisText = {'Str', 'GPe', 'SNr'};
for iRegion = 1:3
    subplot(1, 3, iRegion); hold on;
    scatter(nresponses_per_stim_1([1:end],iRegion), nresponses_per_stim_2([1:end],iRegion),40, 'filled')
    %scatter(responses_per_stim_1(2:3:end,iRegion), responses_per_stim_2(2:3:end,iRegion), 'filled')
    %scatter(responses_per_stim_1(3:3:end,iRegion), responses_per_stim_2(3:3:end,iRegion), 'filled')
    if iRegion ==1
        xlim([0, 0.6])
        ylim([0, 0.6])
    else
        xlim([0, 0.2])
        ylim([0, 0.2])
    end
    
    xlabel(['increase in F.R. after stimulus onset', newline, ...
        '1/2 trials ', thisText{iRegion}])
    ylabel(['increase in F.R. after stimulus onset', newline, ...
        '1/2 trials ', thisText{iRegion}])
    %nan_infs = any(isnan(responses_per_stim_1(1:3:end, :)), 2) | any(isinf(responses_per_stim_2(1:3:end, :)), 2);
    if iRegion ==1
        line([0, 0.6], [0, 0.6])
    else
        line([0, 0.2], [0, 0.2])
    end
    [r, p] = corrcoef(nresponses_per_stim_1(:, 1), nresponses_per_stim_2(:, iRegion), 'rows', 'complete');
    %legend('contralateral', 'central')
prettify_plot;
end


%% scatter: inter-region increase in F.R. 
thisText = {'', 'GPe', 'SNr'};
figure();
for iRegion = 2:3
    subplot(1, 2, iRegion -1)
    hold on;
    %[~, sorted_idx] = sort(responses_per_stim(:, 1));
    scatter(nresponses_per_stim([1:end],1), nresponses_per_stim([1:end],iRegion),40 ,'filled')
    %scatter(nresponses_per_stim([2:3:end],1), responses_per_stim([2:3:end],iRegion), 40, rgb('MidnightBlue'), 'filled')
    %scatter(responses_per_stim(2:3:end, 1), responses_per_stim(2:3:end, iRegion), 'filled') %   0*
    %scatter(responses_per_stim(3:3:end, 1), responses_per_stim(3:3:end, iRegion), 'filled') %  90*
    prettify_plot;
    xlim([0, 0.6])
    ylim([0, 0.6])
    hold on;
    xlabel(['increase in F.R. after stimulus onset', newline, ...
        'striatum'])
    ylabel(['increase in F.R. after stimulus onset', newline, ...
        thisText{iRegion}])

    nan_infs = any(isnan(nresponses_per_stim(:, :)), 2) | any(isinf(nresponses_per_stim(:, :)), 2);
    [r, p] = corrcoef(nresponses_per_stim([1:end], 1), nresponses_per_stim([1:end], iRegion), 'rows', 'complete');
  

    % errorbar(responses_per_stim(second_dataset(:,2)~=90, 1), responses_per_stim(second_dataset(:,2)==-90, iRegion),...
    %     responses_per_stim_std(second_dataset(:,2)~=90, iRegion), "LineStyle", "none", 'Color', 'k')
    % errorbar(responses_per_stim(second_dataset(:,2)~=90, 1), responses_per_stim(second_dataset(:,2)==-90, iRegion),...
    %     responses_per_stim_std(second_dataset(:,2)~=90, 1), 'horizontal', "LineStyle", "none", 'Color', 'k')
    
    line([0, 0.6], [0, 0.6])
    
end


end

