
function [expData, session_data, regions] = cl_loadPerStimulusData(load_type, keep_type, loadVids) %% paramaters
cl_myPaths;

%% get loading info
% recording type exp info
if contains(load_type, 'naive')
    info_table = readtable([csvPath, 'allPassiveRecs.csv'], 'VariableNamingRule', 'modify');
elseif contains(load_type, 'taskGo')
    info_table = readtable([csvPath, 'allTaskGoGoGoRecs.csv'], 'VariableNamingRule', 'modify');
elseif contains(load_type, 'taskNoGo')
    info_table = readtable([csvPath, 'allTaskRecs.csv'], 'VariableNamingRule', 'modify');
end

% regions to load
regions = {'CP', 'GPe', 'SNr'}; % 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
warning off;
regions_id = [672, 1022, 381];
regions_id_old_histo = [574, 611, 823]; % only keep striatum in old ones?, ,611,823];
% get all mice x thisDate x site combinations
%if ismember('OldHisto_', info_table.Properties.VariableNames)
%    use_recs_full = strcmp(info_table.Use_, 'Yes') & ~strcmp(info_table.OldHisto_, 'Yes');
%else
    use_recs_full = strcmp(info_table.Use_, 'Yes');
%end

info_table = sortrows(info_table, 1); % make sure data is sorted by mouse
unique_mice = unique(info_table.Mouse(use_recs_full), 'stable');

% initialize
mouse_thisDate_sites_shank_rec = [];
disp('loading data...')
session_data = struct;
session_data.motion_energy = {};
session_data.trial_type = {};
session_data.session_num = [];

expData = struct;
expData.psth = [];
expData.pvalue = [];
expData.psth_conditions = [];

use_recs = find(use_recs_full);
unitCount = 0;

for iMouse = 1:length(unique_mice)
    thisDate_sites_shank_rec = [info_table.DayNum(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Site(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Shanks(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Rec_Exp(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse)))];

    mouse_thisDate_sites_shank_rec = [mouse_thisDate_sites_shank_rec; repmat(iMouse, size(thisDate_sites_shank_rec, 1), 1), ...
        thisDate_sites_shank_rec];

end


expType = cl_expTypes(load_type);

%% load data

for iRecording = 1:length(use_recs)

    % this recording's info
    animal = unique_mice{mouse_thisDate_sites_shank_rec(iRecording, 1)};
    curr_thisDate = mouse_thisDate_sites_shank_rec(iRecording, 2);
    site = mouse_thisDate_sites_shank_rec(iRecording, 3);
    curr_shank = mouse_thisDate_sites_shank_rec(iRecording, 4);
    curr_rec = mouse_thisDate_sites_shank_rec(iRecording, 5);

    % load probe2ephys info
    probe2ephys_location = cl_cortexlab_filename(animal, [], [], 'probe2ephys');
    if isempty(probe2ephys_location)
        disp('no histo')
        continue;
    end
    load(probe2ephys_location)
    clearvars thisDate_sites_shank_probe
    try
        thisDate_sites_shank_probe = [probe2ephys.day; probe2ephys.site; probe2ephys.shank]';
        old_histo = false;
    catch
        % old histo
        thisDate_sites_shank_probe = [probe2ephys.day; probe2ephys.site; nan(1, size(probe2ephys, 2))]';
        old_histo = true;
    end

    mouse_thisDate_sites_shank_rec(isnan(mouse_thisDate_sites_shank_rec(:, 4)), 4) = 0; %replace nan by 0
    thisDate_sites_shank_probe(isnan(thisDate_sites_shank_probe(:, 3)), 3) = 0; %replace nan by 0

    [~, probe_rec_idx] = ismember(mouse_thisDate_sites_shank_rec(iRecording, 2:4), thisDate_sites_shank_probe, 'rows');

    % get probe_ccf
    probe_ccf_location = cl_cortexlab_filename(animal, [], [], 'histo');
    load(probe_ccf_location)
    new_units = [];
    units_to_keep = [];

    % get all experiments
    experiments = cl_find_experiments(animal, '', true);
    experiments = experiments([experiments.ephys]);

    % are any of the ROIs present
    this_probe = probe_rec_idx;
    if this_probe > 0
        for iRegion = 1:size(regions, 2)
            if size(probe_ccf, 1) < this_probe
                % Get field names
                fields = fieldnames(probe_ccf);
                
                % Iterate over each field and add the value
                for i = 1:length(fields)
                    fieldName = fields{i}; % Get field name
                    probe_ccf(this_probe).(fieldName) = 0; % Add value
                end
            end
                if old_histo
                    this_region = regions_id_old_histo(iRegion);
                    this_region_idx = find(probe_ccf(this_probe).trajectory_areas == this_region);
                else

                    this_region = regions_id(iRegion);
                    this_region_idx = find(probe_ccf(this_probe).trajectory_areas == this_region);

                end
           
            if ~isempty(this_region_idx)
                try
                    this_region_start(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(1));
                    this_region_stop(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(end));
                catch
                    this_region_start(iRegion) = 0;
                    this_region_stop(iRegion) = 0;
                end
            else
                this_region_start(iRegion) = 0;
                this_region_stop(iRegion) = 0;
            end
        end

        % load units in relevant depths
        these_regions_present = find(this_region_stop);
        if ~isempty(these_regions_present)
            these_exps = str2num(info_table.Exps{use_recs(iRecording)});

            if isnan(curr_rec)
                recording = [];
            else
                recording = curr_rec;
            end

            goodExp_max_trials = cl_get_max_good_experiment(experiments, recording, site, these_exps, curr_thisDate, animal, expType{keep_type});
            if isempty(goodExp_max_trials)
                continue;
            end
            for iExperiment = goodExp_max_trials
                experiment = these_exps(iExperiment);
                experiments = cl_find_experiments(animal, '', true);
                experiments = experiments([experiments.ephys]);
                thisDate = experiments(curr_thisDate).thisDate;

                verbose = false; % display load progress and some info figures
                if loadVids
                    load_parts.cam = true;
                end
                load_parts.imaging = false;
                load_parts.ephys = true;
                loadClusters = 0;
                try
                    cl_load_experiment;
                catch
                    continue;
                end


                %subselect shank
                if curr_shank > 0
                    shank_xdepth = [250 * (curr_shank - 1), 250 * (curr_shank - 1) + 32];
                    shank_units = find(template_xdepths >= shank_xdepth(1) & template_xdepths <= shank_xdepth(2));
                else
                    shank_units = 1:size(template_xdepths, 1)';
                end
                % get units we want to keep + store each units location (region + coordinates)
                units_to_keep = [];
                units_to_keep_area = [];
                units_to_keep_coords = [];
                unit_side = [];

                for iRegion = 1:length(these_regions_present)
                    new_units = find(template_depths(shank_units) >= this_region_start(these_regions_present(iRegion)) & ...
                        template_depths(shank_units) <= this_region_stop(these_regions_present(iRegion)));
                    units_to_keep = [units_to_keep; new_units];
                    units_to_keep_area = [units_to_keep_area; ones(size(new_units, 1), 1) .* these_regions_present(iRegion)];
                    %units_to_keep_coords

                    unit_closest_depth = arrayfun(@(x) ... %closest depth
                        find(probe_ccf(this_probe).probe_depths >= template_depths(shank_units(new_units(x))), 1, 'first'), 1:length(new_units));
                    if old_histo
                        units_to_keep_coords = [squeeze(units_to_keep_coords); ...
                            squeeze(probe_ccf(this_probe).trajectory_coords(unit_closest_depth, :))]./2.5;
                        % try
                        %     units_coords_regularized = [probe_ccf(this_probe).trajectory_coords(unit_closest_depth, 1), ...
                        %         probe_ccf(this_probe).trajectory_coords(unit_closest_depth, 3), ...
                        %         probe_ccf(this_probe).trajectory_coords(unit_closest_depth, 2)]./2.5;
                        %     units_to_keep_coords = [squeeze(units_to_keep_coords); ...
                        %     squeeze(units_coords_regularized)]; %2.5 scaling factor bewteen allen and brainreg
                        % catch
                        %     keyboard;
                        % end
                        
                    else
                        units_to_keep_coords = [units_to_keep_coords; ...
                            probe_ccf(this_probe).trajectory_coords(unit_closest_depth, :)];
                    end
                    % AP, DV, ML
                    bregma = [540, 0, 570];
                    unit_side = [unit_side; ...
                        probe_ccf(this_probe).trajectory_coords(unit_closest_depth, 3) - bregma(1) / 2.5 < 0]; %-1 for left, 1 for right
                end
                protocol = '';
                rerunEP = 0;
                plotGUI = 0;
                runQM = 1;
                rerunQM = 0;
                region = '';
                runEP = 1;
                clearvars unitType
                try
                    [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, thisDate, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
                    %bc_qualityMetricsPipeline_JF(animal, thisDate, site, recording, experiment_num, protocol, rerunQM, plotGUI, runQM)
                catch
                    rerunQM = 1;
                    [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, thisDate, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
                end

                ephysProperties = bc_ephysPropertiesPipeline_JF(animal, thisDate, site, recording, 1, rerunEP, runEP, region);


                expData.psth_conditions_all{iRecording, keep_type} = unique(trial_conditions, 'rows');

                % get stim response
                unique_templates = unique(spike_templates);
                if contains(load_type, 'task')
                    if keep_type == 1
                        expData.psth_conditions{keep_type} = unique(trial_conditions, 'rows');

                    else
                        expData.psth_conditions{keep_type} = unique(trial_conditions(ismember(trial_conditions(:, 1), ...
                            [1:13]) & ismember(trial_conditions(:, 2), [-90, 0, 90]), :), 'rows');
                    end
                else
                    %if keep_type == 1%nat images
                    expData.psth_conditions{keep_type} = unique(trial_conditions, 'rows');
                    expData.psth_conditions_all{iRecording, keep_type} = unique(trial_conditions, 'rows');
                    if keep_type == 5
                        if size(trial_conditions, 2) == 2

                            these_guys = ismember(trial_conditions(:, 1), ...
                                [1:13]) & ismember(trial_conditions(:, 2), [-90, 0, 90]);
                            disp(unique(trial_conditions, 'rows'))

                            %if size(unique(trial_conditions(these_guys,:), 'rows'),1) ==3
                            %passive_data_per_cond.psth_conditions{keep_type} = unique(trial_conditions(these_guys,:), 'rows');
                            %else
                            %    disp('no stim screen location saved')
                            %    continue;
                            %end
                        else
                            continue;
                        end
                    else
                        %passive_data_per_cond.psth_conditions_all{iRecording, keep_type} = unique(trial_conditions, 'rows');
                    end


                end

                % for half trials, identify max.
                % then average on that half.
                expData.trial_types{keep_type, iRecording} = trial_conditions;
                if exist('no_move_trials', 'var')
                    expData.no_move_trials{keep_type, iRecording} = no_move_trials;
                end
                raster_window = [-0.25, 0.6];

                if loadVids
                    try
                        psth_bin_size = 0.01;
                        camera_location = facecam_fn;
                        camera_t = facecam_t;
                        framerate = 30;
                        surround_t = raster_window;
                        use_align = stimOn_times;
                        surround_frames = diff(raster_window) * framerate;

                        use_t_plotting = '';
                        plot_me = false;
                        box_crop = [];

                        [motion_energy, motion_energy_std, cam_align_diff_avg, n_trials, motion_energy_frame_full] = cl_getMotionIndex(camera_location, ...
                            camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

                        session_data.motion_energy{iRecording} = motion_energy_frame_full;
                        session_data.trial_type{iRecording} = trial_conditions;
                        session_data.session_num = iRecording;
                    catch
                        session_data.motion_energy{iRecording} = NaN;
                        session_data.trial_type{iRecording} = NaN;
                        session_data.session_num = iRecording;
                    end
                end

                psth_bin_size = 0.001;
                for iUnit = 1:size(units_to_keep, 1)
                    if exist('no_move_trials', 'var') % only in passive protocols
                        keep_trials = no_move_trials;
                    else
                        keep_trials = ones(n_trials(end), 1);
                    end


                    align_times = stimOn_times(keep_trials);
                    [~, trial_cond_idx] = ismember(trial_conditions(keep_trials, :), expData.psth_conditions{keep_type}, 'rows');


                    [~, curr_raster, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times, []);

                    % p value test for *each* condition
                    for iCond = 1:size(expData.psth_conditions{keep_type}, 1)
                        if ~isempty(keep_trials)
                            [~, trial_cond_idx_single] = ismember(trial_conditions(keep_trials, :), expData.psth_conditions{keep_type}(iCond, :), 'rows');
                        else
                            [~, trial_cond_idx_single] = ismember(trial_conditions, expData.psth_conditions{keep_type}(iCond, :), 'rows');

                        end
                        pvals_per_cond(iCond) = signrank(nanmean(curr_raster(logical(trial_cond_idx_single), 1:150), 2), ...
                            nanmean(curr_raster(logical(trial_cond_idx_single), 250:400), 2));
                        expData.pvalue{keep_type}(iUnit + unitCount, iCond) = pvals_per_cond(iCond);


                        pre_activity = nanmean(curr_raster(logical(trial_cond_idx_single), 1:150), 2);
                        post_activity = nanmean(curr_raster(logical(trial_cond_idx_single), 301:450), 2);
                        rand_activity = nanmean(curr_raster(logical(trial_cond_idx_single), 151:300), 2);
                        all_activity = [pre_activity; post_activity];
                        shuffling_idx = randi(size(all_activity, 1), 1000);
                        shuffled_diffs = arrayfun(@(x) nanmean(all_activity(shuffling_idx(1:500, x)) ...
                            -all_activity(shuffling_idx(501:1000, x))), 1:1000);
                        real_diff = nanmean(post_activity-pre_activity);
                        fake_diff = nanmean(rand_activity-pre_activity);
                        pctile_2 = prctile(shuffled_diffs, 95);
                        pctile_1 = prctile(shuffled_diffs, 5);

                        expData.pvalue_shuffled_005_per_cond{keep_type}(iUnit + unitCount, iCond) = ...
                            real_diff >= pctile_2 || real_diff <= pctile_1;

                        expData.pvalue_shuffled_0025_per_cond{keep_type}(iUnit + unitCount, iCond) = ...
                            real_diff >= prctile(shuffled_diffs, 97.5) || real_diff <= prctile(shuffled_diffs, 2.5);


                         expData.pvalue_shuffled_005_per_cond_fake{keep_type}(iUnit + unitCount) = ...
                            fake_diff >= pctile_2 || fake_diff <= pctile_1;

                        expData.pvalue_shuffled_0025_per_cond_fake{keep_type}(iUnit + unitCount) = ...
                            fake_diff >= prctile(shuffled_diffs, 97.5) || fake_diff <= prctile(shuffled_diffs, 2.5);
                        
                        clearvars curr_raster_smoothed_fr
                        curr_raster_smoothed_fr = smoothdata(curr_raster(logical(trial_cond_idx_single),:), 2, 'gaussian', [0 50]).*(1/psth_bin_size);
                 

                        pre_activity_max = nanmax(curr_raster_smoothed_fr(:, 1:150),[],2);
                        post_activity_max = nanmax(curr_raster_smoothed_fr(:,301:450),[],2);
                        rand_activity_max = nanmax(curr_raster_smoothed_fr(:, 151:300),[],2);
                        all_activity_max = [pre_activity_max; post_activity_max];
                        shuffling_idx = randi(size(all_activity_max, 1), 1000);
                        shuffled_diffs = arrayfun(@(x) nanmean(all_activity_max(shuffling_idx(1:500, x)) ...
                            -all_activity_max(shuffling_idx(501:1000, x))), 1:1000);
                        real_diff = nanmean(post_activity_max-pre_activity_max);
                        fake_diff = nanmean(rand_activity_max-pre_activity_max);
                        pctile_2 = prctile(shuffled_diffs, 95);
                        pctile_1 = prctile(shuffled_diffs, 5);

                        expData.pvalue_shuffled_005_per_cond_peak{keep_type}(iUnit + unitCount, iCond) = ...
                            real_diff >= pctile_2 ;

                        expData.pvalue_shuffled_0025_per_cond_peak{keep_type}(iUnit + unitCount, iCond) = ...
                            real_diff >= prctile(shuffled_diffs, 97.5) ;

                        expData.pvalue_shuffled_005_per_cond_peak_fake{keep_type}(iUnit + unitCount) = ...
                            fake_diff >= pctile_2;

                        expData.pvalue_shuffled_0025_per_cond_peak_fake{keep_type}(iUnit + unitCount) = ...
                            fake_diff >= prctile(shuffled_diffs, 97.5) ;


                        pre_activity_min = nanmin(curr_raster_smoothed_fr(:, 1:150),[],2);
                        post_activity_min = nanmin(curr_raster_smoothed_fr(:,301:450),[],2);
                        rand_activity_min = nanmin(curr_raster_smoothed_fr(:, 151:300),[],2);
                        all_activity_min = [pre_activity_min; post_activity_min];
                        shuffling_idx = randi(size(all_activity_min, 1), 1000);
                        shuffled_diffs = arrayfun(@(x) nanmean(all_activity_min(shuffling_idx(1:500, x)) ...
                            -all_activity_min(shuffling_idx(501:1000, x))), 1:1000);
                        real_diff = nanmean(post_activity_min-pre_activity_min);
                        fake_diff = nanmean(rand_activity_min-pre_activity_min);
                        pctile_2 = prctile(shuffled_diffs, 95);
                        pctile_1 = prctile(shuffled_diffs, 5);

                        expData.pvalue_shuffled_005_per_cond_trough{keep_type}(iUnit + unitCount, iCond) = ...
                            real_diff <= pctile_1;

                        expData.pvalue_shuffled_0025_per_cond_trough{keep_type}(iUnit + unitCount, iCond) = ...
                             real_diff <= prctile(shuffled_diffs, 2.5);

                        expData.pvalue_shuffled_005_per_cond_trough_fake{keep_type}(iUnit + unitCount, iCond) = ...
                            fake_diff <= pctile_1;

                        expData.pvalue_shuffled_0025_per_cond_trough_fake{keep_type}(iUnit + unitCount, iCond) = ...
                            fake_diff <= prctile(shuffled_diffs, 2.5);


                    end


                    % across all conds shuffle test
                    pre_activity = nanmean(curr_raster(:, 1:150), 2);
                    post_activity = nanmean(curr_raster(:, 301:450), 2);
                    rand_activity = nanmean(curr_raster(:, 151:300), 2);
                    all_activity = [pre_activity; post_activity];
                    shuffling_idx = randi(size(all_activity, 1), 1000);
                    shuffled_diffs = arrayfun(@(x) nanmean(all_activity(shuffling_idx(1:500, x)) ...
                        -all_activity(shuffling_idx(501:1000, x))), 1:1000);
                    real_diff = nanmean(post_activity-pre_activity);
                    fake_diff = nanmean(rand_activity-pre_activity);
                    pctile_2 = prctile(shuffled_diffs, 95);
                    pctile_1 = prctile(shuffled_diffs, 5);

                    expData.pvalue_shuffled_005{keep_type}(iUnit + unitCount) = ...
                        real_diff >= pctile_2 || real_diff <= pctile_1;

                    expData.pvalue_shuffled_0025{keep_type}(iUnit + unitCount) = ...
                        real_diff >= prctile(shuffled_diffs, 97.5) || real_diff <= prctile(shuffled_diffs, 2.5);

                    expData.pvalue_shuffled_005_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff >= pctile_2 || fake_diff <= pctile_1;

                    expData.pvalue_shuffled_0025_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff >= prctile(shuffled_diffs, 97.5) || fake_diff <= prctile(shuffled_diffs, 2.5);


                    % smooth raster, then take f.r.
                    curr_raster_smoothed_fr = smoothdata(curr_raster, 2, 'gaussian', [0 50]).*(1/psth_bin_size);
                    pre_activity = nanmax(curr_raster_smoothed_fr(:, 1:150), [], 2);
                    post_activity = nanmax(curr_raster_smoothed_fr(:, 301:450),[],  2);
                    rand_activity = nanmax(curr_raster_smoothed_fr(:, 151:300),[],  2);
                    all_activity = [pre_activity; post_activity];
                    shuffling_idx = randi(size(all_activity, 1), 1000);
                    shuffled_diffs = arrayfun(@(x) nanmean(all_activity(shuffling_idx(1:500, x)) ...
                        -all_activity(shuffling_idx(501:1000, x))), 1:1000);
                    real_diff = nanmean(post_activity-pre_activity);
                    fake_diff = nanmean(rand_activity-pre_activity);
                    pctile_2 = prctile(shuffled_diffs, 95);

                    expData.pvalue_shuffled_005_peak{keep_type}(iUnit + unitCount) = ...
                        real_diff >= pctile_2;

                    expData.pvalue_shuffled_0025_peak{keep_type}(iUnit + unitCount) = ...
                        real_diff >= prctile(shuffled_diffs, 97.5);

                    expData.pvalue_shuffled_005_peak_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff >= pctile_2;

                    expData.pvalue_shuffled_0025_peak_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff >= prctile(shuffled_diffs, 97.5);


                    pre_activity = nanmin(curr_raster_smoothed_fr (:, 1:150), [], 2);
                    post_activity = nanmin(curr_raster_smoothed_fr (:, 301:450),[],  2);
                    rand_activity = nanmin(curr_raster_smoothed_fr (:, 151:300), [], 2);
                    all_activity = [pre_activity; post_activity];
                    shuffling_idx = randi(size(all_activity, 1), 1000);
                    shuffled_diffs = arrayfun(@(x) nanmean(all_activity(shuffling_idx(1:500, x)) ...
                        -all_activity(shuffling_idx(501:1000, x))), 1:1000);
                    real_diff = nanmean(post_activity-pre_activity);
                    fake_diff = nanmean(rand_activity-pre_activity);
                    pctile_1 = prctile(shuffled_diffs, 5);

                    expData.pvalue_shuffled_005_trough{keep_type}(iUnit + unitCount) = ...
                        real_diff <= pctile_1;

                    expData.pvalue_shuffled_0025_trough{keep_type}(iUnit + unitCount) = ...
                        real_diff <= prctile(shuffled_diffs, 2.5);

                     expData.pvalue_shuffled_005_trough_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff <= pctile_1;

                    expData.pvalue_shuffled_0025_trough_fake{keep_type}(iUnit + unitCount) = ...
                        fake_diff <= prctile(shuffled_diffs, 2.5);

                    
                    
                    raster_window = [-0.2, 0.6];
                    % psth half trials 1
                    [curr_psth, ~, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(1:2:end), trial_cond_idx(1:2:end));
                    expData.av_psth_1{keep_type, iRecording}(iUnit, :, :) = curr_psth;


                    expData.psth{keep_type}(iUnit + unitCount, 1, 1:size(curr_psth, 1), :) = curr_psth;

                    % psth half trials 2
                    [curr_psth, ~, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(2:2:end), trial_cond_idx(2:2:end));

                    expData.psth{keep_type}(iUnit + unitCount, 2, 1:size(curr_psth, 1), :) = curr_psth;
                    expData.av_psth_2{keep_type, iRecording}(iUnit, :, :) = curr_psth;


                    % psth all trials
                    [curr_psth, curr_raster, t_det, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(1:1:end), trial_cond_idx(1:1:end));
                    expData.psth{keep_type}(iUnit + unitCount, 3, 1:size(curr_psth, 1), :) = curr_psth;
                    startIdx = find(t_det >= 0.05, 1, 'first');
                    stopIdx = find(t_det >= 0.2, 1, 'first');
                    expData.av_per_trial{keep_type, iRecording}(iUnit, :) = nanmean(curr_raster(:, startIdx:stopIdx), 2);
                    startIdx = find(t_det >= -0.2, 1, 'first');
                    stopIdx = find(t_det >= -0.05, 1, 'first');
                    expData.av_per_trial_base{keep_type, iRecording}(iUnit, :) = nanmean(curr_raster(:, startIdx:stopIdx), 2);


                end

            end

            % save data in structure
            if ~isempty(units_to_keep)
                expData.nz_trial_types{keep_type, iRecording} = trial_conditions;
                expData.nz_no_move_trials{keep_type, iRecording} = no_move_trials;
                expData.animal_thisDate_site_shank(unitCount+1:unitCount+size(units_to_keep, 1), :) = ...
                    repmat([mouse_thisDate_sites_shank_rec(iRecording, 1), curr_thisDate, site, curr_shank], size(units_to_keep, 1), 1);
                expData.unit_area(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_area;
                expData.unit_coords(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_coords;
                expData.unit_side(unitCount+1:unitCount+size(units_to_keep, 1)) = unit_side;
                expData.t = t;
                %                passive_data_per_cond.t_det = t_det;
                expData.unitNum((unitCount + 1:unitCount + size(units_to_keep, 1))) = shank_units(units_to_keep);
                %try
                try
                    expData.unitType((unitCount + 1:unitCount + size(units_to_keep, 1))) = unitType(units_to_keep);
                catch
                    rerunQM = 1;
                    [unitType, ~] = bc_qualityMetricsPipeline_JF(animal, thisDate, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
                    expData.unitType((unitCount + 1:unitCount + size(units_to_keep, 1))) = unitType(units_to_keep);
                end
                
                if isstruct(ephysProperties)
                    ephysProperties = struct2table(ephysProperties);
                end
                if ismember('propLongISI', ephysProperties.Properties.VariableNames)
                    if size(ephysProperties.propLongISI, 1) == size(unitType, 1)
                        try
                            expData.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression(units_to_keep);
                            expData.templateDuration((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.templateDuration(units_to_keep);
                            expData.propLongISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                            expData.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.spike_rateSimple(units_to_keep);
                        catch
                            expData.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression_ms(units_to_keep);
                            expData.templateDuration((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.waveformDuration_peakTrough_us(units_to_keep);
                            expData.propLongISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                            expData.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.mean_firingRate(units_to_keep);
                        end
                    else

                        rerunEP = 1;
                        ephysProperties = bc_ephysPropertiesPipeline_JF(animal, thisDate, site, recording, 1, rerunEP, runEP, region);
                        expData.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression_ms(units_to_keep);
                        expData.templateDuration((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.waveformDuration_peakTrough_us(units_to_keep);
                        expData.propLongISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                        expData.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.mean_firingRate(units_to_keep);

                    end
                else

                    rerunEP = 1;
                    ephysProperties = bc_ephysPropertiesPipeline_JF(animal, thisDate, site, recording, 1, rerunEP, runEP, region);
                    expData.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression_ms(units_to_keep);
                    expData.templateDuration((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.waveformDuration_peakTrough_us(units_to_keep);
                    expData.propLongISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                    expData.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.mean_firingRate(units_to_keep);

                end

                %passive_data_per_cond
                unitCount = unitCount + size(units_to_keep, 1);

            end

            % clear variables
            disp(['   ', num2str(iRecording), '/', num2str(length(use_recs))])

            keep regions_id_old_histo expType session_data mouse_thisDate_sites_shank_rec unique_mice expData use_recs info_table regions regions_id unitCount load_type keep_type loadVids
        end
    end
    %catch
    %end
end

end