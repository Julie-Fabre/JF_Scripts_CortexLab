%keepVis = 0;
%passive = 0;
%goNogo = 1;
keep passive goNogo keepVis

if passive
    task_data_here = load('/home/julie/Dropbox/MATLAB/task_data_passive.mat');
    idx = 5; %2: task, 4/5 : passive

elseif goNogo
    task_data_here = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
    idx = 2; %2: task, 4/5 : passive
else
    task_data_here = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
    idx = 2; %2: task, 4/5 : passive

end

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

for iRegion = 1:size(plot_regions, 2)
    if keepVis
        these_units = task_data_here.unit_area == iRegion & ...
            (task_data_here.unitType' == 1 | task_data_here.unitType' == 2) & ...
            task_data_here.pvalue_shuffled_005{1, idx}' == 1;
        these_units_vis = task_data_here.unit_area == iRegion & ...
            (task_data_here.unitType' == 1 | task_data_here.unitType' == 2) & ...
            task_data_here.pvalue_shuffled_005{1, idx}' == 1;
    else
        these_units = task_data_here.unit_area == iRegion & ...
            (task_data_here.unitType' == 1 | task_data_here.unitType' == 2); % & ...
        %task_data_here.pvalue_shuffled_005{1,idx}' == 1;
          these_units_vis = task_data_here.unit_area == iRegion & ...
            (task_data_here.unitType' == 1 | task_data_here.unitType' == 2) & ...
            task_data_here.pvalue_shuffled_005{1, idx}' == 1;
    end
    unitCount = 0;
    for iSession = 1:size(nonZero_idx, 2) - 1
        thisSession = nonZero_idx(iSession+1) - 1;
        trials_no_move = task_data_here.trial_types{idx, thisSession}(task_data_here.no_move_trials{idx, thisSession}, :);
        %these_units_session = these_units(cumsum_cells(iSession):cumsum_cells(iSession+1));
        these_units_session = these_units(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);
        these_units_session_vis = these_units_vis(session_cumCells(thisSession):session_cumCells(thisSession+1)-1);

        if ~isempty(trials_no_move)
             activity_per_trial_neuron = task_data_here.av_per_trial{idx, thisSession}(these_units_session, :)';
          
            if ~isempty(activity_per_trial_neuron) && ~isempty(find(any(activity_per_trial_neuron > 0)))
                
                increaseFR_session_fraction(iRegion, iSession) = sum(these_units_session_vis) ./ length(these_units_session);
            else
                increaseFR_session_fraction(iRegion, iSession) = NaN;
            end

        else
           increaseFR_session_fraction(iRegion, iSession) = NaN;
        end

    end

end