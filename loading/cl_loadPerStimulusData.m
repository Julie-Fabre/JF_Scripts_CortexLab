function passive_data_per_cond = cl_loadPerStimulusData(load_type, keep_type) %% paramaters
cl_myPaths;
    if contains(load_type, 'passive')
    info_table = readtable([csvPath 'allPassiveRecs.csv'], 'VariableNamingRule', 'modify');
elseif contains(load_type, 'taskGo')
    info_table = readtable([csvPath 'old' filesep 'allTaskGoGoGoRecs_prev.csv'], 'VariableNamingRule', 'modify');
else
    info_table =readtable([csvPath 'allTaskRecs.csv'], 'VariableNamingRule', 'modify');
    end

load_passive_ = true;
%info_table = readtable('/home/julie/Dropbox/PhD_summary/allPassiveRecs.csv');
regions = {'CP', 'GPe', 'SNr'};% 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
%regions = {'CP'};
warning off;
%get region nums (allen format)

if ~exist('st', 'var')
   % [tv, av, st, bregma] = bd_loadAllenAtlas(atlasBrainRegLocation);
    st =[];
end
regions_id = [672, 1022, 381];
% for iRegion = 1:size(regions, 2)
%     regions_id(iRegion) = st.id(strcmp(st.acronym, regions(iRegion)));
% end

% get all mice x day x site combinations
use_recs_full = strcmp(info_table.Use_, 'Yes');

% good units (bombcell) (+ cell type?) . location. psth per image type.

%% load all mouse rec locations
info_table = sortrows(info_table, 1); % make sure data is sorted by mouse
unique_mice = unique(info_table.Mouse(use_recs_full));

% initialize
mouse_day_sites_shank_rec = [];
disp('loading data...')
passive_data_per_cond = struct;
passive_data_per_cond.psth = [];
passive_data_per_cond.pvalue = [];
passive_data_per_cond.psth_conditions = [];

use_recs = find(use_recs_full);
unitCount = 0;

for iMouse = 1:length(unique_mice)
    day_sites_shank_rec = [info_table.DayNum(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Site(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Shanks(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Rec_Exp(use_recs_full & strcmp(info_table.Mouse, unique_mice(iMouse)))];

    mouse_day_sites_shank_rec = [mouse_day_sites_shank_rec; repmat(iMouse, length(day_sites_shank_rec), 1), ...
        day_sites_shank_rec];

end

%% load data

for iRecording = 1:length(use_recs)%61:length(use_recs)
    %try
    %try
    % curr variables

    animal = unique_mice{mouse_day_sites_shank_rec(iRecording, 1)};
    curr_day = mouse_day_sites_shank_rec(iRecording, 2);
    site = mouse_day_sites_shank_rec(iRecording, 3);
    curr_shank = mouse_day_sites_shank_rec(iRecording, 4);
    curr_rec = mouse_day_sites_shank_rec(iRecording, 5);

    % load probe_ccf and probe2ephys
    probe2ephys_location = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys');
    if isempty(probe2ephys_location)
        disp('no histo')
        continue;
    end
    load(probe2ephys_location)
    day_sites_shank_probe = [probe2ephys.day; probe2ephys.site; probe2ephys.shank]';

    mouse_day_sites_shank_rec(isnan(mouse_day_sites_shank_rec(:, 4)), 4) = 0; %replace nan by 0
    day_sites_shank_probe(isnan(day_sites_shank_probe(:, 3)), 3) = 0; %replace nan by 0

    [probe_rec, probe_rec_idx] = ismember(mouse_day_sites_shank_rec(iRecording, 2:4), day_sites_shank_probe, 'rows');

    probe_ccf_location = AP_cortexlab_filenameJF(animal, [], [], 'histo');
    load(probe_ccf_location)
    new_units =[];
    units_to_keep =[];
    % is any region present
    this_probe = probe_rec_idx;
    if this_probe > 0
        for iRegion = 1:size(regions, 2)
            this_region = regions_id(iRegion);
            this_region_idx = find(probe_ccf(this_probe).trajectory_areas == this_region);
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
            % in order: gratings, location, (exp type 1,2,3,4)
            % conditions = trial_conditions
            if contains(load_type, 'passive')
                load_me = 1:size(these_exps, 2);
            else
                load_me = size(these_exps, 2);
            end
            % find max experiment 
            n_trials = [];
            exp_n_trials = [];
            for iExperiment = load_me
                experiment = these_exps(iExperiment); %for now, just load 1 experiment
                experiments = AP_find_experimentsJF(animal, '', true);
                experiments = experiments([experiments.ephys]);
                day = experiments(curr_day).day;
                if isnan(curr_rec)
                    recording = [];
                else
                    recording = curr_rec;
                end


                verbose = false; % display load progress and some info figures
                load_parts.cam = false;
                load_parts.imaging = false;
                load_parts.ephys = true;
                loadClusters = 0;
                try
                    JF_load_experiment;
                        
                catch
                    disp('error loading experiment')
                    continue;
                end
                if contains(load_type, 'passive')
                    if contains(expDef, 'JF_GratingPassiveVarITI')
                        if keep_type ~= 1
                        continue;
                        end
                        experiment_type = 1;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_locations')
                        if keep_type ~= 2
                        continue;
                        end
                        experiment_type = 2;
                    elseif contains(expDef, 'JF_natural_imagesFit')
                        if keep_type ~= 3
                        continue;
                        end
                        experiment_type = 3;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_choiceworldStimuli_onlyTask')
                        if keep_type ~= 4
                        continue;
                        end
                        experiment_type = 4;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_choiceworldStimuli')
                        if keep_type ~= 5
                        continue;
                        end
                        experiment_type = 5;
                        keep_trial = no_move_trials;
                        
                    elseif contains(expDef, 'JF_natural_images')
                        if keep_type ~= 6
                        continue;
                        end
                        experiment_type = 6;
                        keep_trial = no_move_trials;
                    end
                else
                    if contains(expDef, 'noGo_')
                        if keep_type ~= 1
                            continue;
                        end
                        experiment_type = 1;
                        keep_trial = [];
                    elseif contains(expDef, 'JF_choiceworldStimuli')
                        if keep_type ~= 2
                            continue;
                        end
                        experiment_type = 2;
                        keep_trial = no_move_trials;
                    end
                end
                n_trials(iExperiment) = length(keep_trial);
                exp_n_trials(iExperiment) = experiment;

            end
            
            load_me_max = exp_n_trials(find(n_trials==max(n_trials)));
            for iExperiment = load_me_max %1:size(these_exps,2)
                experiment =  load_me_max;%these_exps(iExperiment); %for now, just load 1 experiment
                experiments = AP_find_experimentsJF(animal, '', true);
                experiments = experiments([experiments.ephys]);
                day = experiments(curr_day).day;
                if isnan(curr_rec)
                    recording = [];
                else
                    recording = curr_rec;
                end


                verbose = false; % display load progress and some info figures
                load_parts.cam = false;
                load_parts.imaging = false;
                load_parts.ephys = true;
                loadClusters = 0;
                try
                    JF_load_experiment;
                        
                catch
                    disp('error loading experiment')
                    continue;
                end
                if contains(load_type, 'passive')
                    if contains(expDef, 'JF_GratingPassiveVarITI')
                        if keep_type ~= 1
                        continue;
                        end
                        experiment_type = 1;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_locations')
                        if keep_type ~= 2
                        continue;
                        end
                        experiment_type = 2;
                    elseif contains(expDef, 'JF_natural_imagesFit')
                        if keep_type ~= 3
                        continue;
                        end
                        experiment_type = 3;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_choiceworldStimuli_onlyTask')
                        if keep_type ~= 4
                        continue;
                        end
                        experiment_type = 4;
                        keep_trial = no_move_trials;
                    elseif contains(expDef, 'JF_choiceworldStimuli')
                        if keep_type ~= 5
                        continue;
                        end
                        experiment_type = 5;
                        keep_trial = no_move_trials;
                        
                    elseif contains(expDef, 'JF_natural_images')
                        if keep_type ~= 6
                        continue;
                        end
                        experiment_type = 6;
                        keep_trial = no_move_trials;
                    end
                else
                    if contains(expDef, 'noGo_')
                        if keep_type ~= 1
                            continue;
                        end
                        experiment_type = 1;
                        keep_trial = [];
                    elseif contains(expDef, 'JF_choiceworldStimuli')
                        if keep_type ~= 2
                            continue;
                        end
                        experiment_type = 2;
                        keep_trial = no_move_trials;
                    end
                end
                loadClusters = 0;
                % "real" template depths
                %         metaFile = strrep(ephysAP_path, '.cbin', '.meta');
                %         [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);
                %         if contains(channelMapIMRO, 'hStripe')
                %             template_depths = template_depths - 2880 + 750;
                %         end

                %subselect shank
                if ismember(experiment_type, keep_type)
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
                    
                    for iRegion = 1:length(these_regions_present)
                        new_units = find(template_depths(shank_units) >= this_region_start(these_regions_present(iRegion)) & ...
                            template_depths(shank_units) <= this_region_stop(these_regions_present(iRegion)));
                        units_to_keep = [units_to_keep; new_units];
                        units_to_keep_area = [units_to_keep_area; ones(size(new_units, 1), 1) .* these_regions_present(iRegion)];

                        unit_closest_depth = arrayfun(@(x) ... %closest depth
                            find(probe_ccf(this_probe).probe_depths >= template_depths(shank_units(new_units(x))), 1, 'first'), 1:length(new_units));
                        units_to_keep_coords = [units_to_keep_coords; ...
                            probe_ccf(this_probe).trajectory_coords(unit_closest_depth, :)];
                    end
                    protocol = '';
                    rerunEP = 0;
                    plotGUI = 0;
                    runQM = 1;
                    rerunQM = 0;
                    region = '';
                    runEP = 1;
                    try
                        [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, day, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
                    catch
                        rerunQM = 1;
                        [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, day, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
                    end
                    ephysProperties = bc_ephysPropertiesPipeline_JF(animal, day, site, recording, 1, protocol, rerunEP, runEP, region);

                end


                % get stim response
                unique_templates = unique(spike_templates);
                if contains(load_type, 'task')
                    if experiment_type == 1
                        passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions, 'rows');
                    
                   else
                        passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions(ismember(trial_conditions(:, 1), ...
                            [1:13]) & ismember(trial_conditions(:, 2), [-90, 0,90]), :), 'rows');
                    end
                else
                    %if experiment_type == 1%nat images
                    passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions, 'rows');
                    if experiment_type==5
                        if size(trial_conditions,2) == 2
                     
                            these_guys = ismember(trial_conditions(:, 1), ...
                                [4,6,12]) & ismember(trial_conditions(:, 2), [-90]);
                            disp(unique(trial_conditions(these_guys,:), 'rows'))
                            if size(unique(trial_conditions(these_guys,:), 'rows'),1) ==3
                                passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions(these_guys,:), 'rows');
                            else
                                disp('no stim screen location saved')
                                continue;
                            end
                        else
                            continue;
                        end
                    else
                        passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions, 'rows');
                    end
                        

                end

                % for half trials, identify max.
                % then average on that half.
                passive_data_per_cond.trial_types{experiment_type, iRecording} = trial_conditions;
                passive_data_per_cond.no_move_trials{experiment_type, iRecording} = no_move_trials;
                for iUnit = 1:size(units_to_keep, 1)
                    raster_window = [-0.5, 1];
                    psth_bin_size = 0.01;
                    if ~isempty(keep_trial)
                        align_times = stimOn_times(keep_trial);
                        [~, trial_cond_idx] = ismember(trial_conditions(keep_trial,:), passive_data_per_cond.psth_conditions{experiment_type}, 'rows');
                    
                    else
                        align_times = stimOn_times;
                        [~, trial_cond_idx] = ismember(trial_conditions, passive_data_per_cond.psth_conditions{experiment_type}, 'rows');
                    
                    end

                    %stats = grpstats(passive_data_per_cond.psth_conditions{experiment_type}, trial_conditions);
                    %responsive cell


                    [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times, []);
                    % p value test for *each* condition 
                    for iCond = 1:size(passive_data_per_cond.psth_conditions{experiment_type},1)
                         if ~isempty(keep_trial)
                             [~, trial_cond_idx_single] = ismember(trial_conditions(keep_trial,:), passive_data_per_cond.psth_conditions{experiment_type}(iCond,:), 'rows');
                         else
                             [~, trial_cond_idx_single] = ismember(trial_conditions, passive_data_per_cond.psth_conditions{experiment_type}(iCond,:), 'rows');
                        
                         end
                        pvals_per_cond(iCond) = signrank(nanmean(curr_raster(logical(trial_cond_idx_single), 40:50), 2),...
                            nanmean(curr_raster(logical(trial_cond_idx_single), 55:65), 2));
                         passive_data_per_cond.pvalue{experiment_type}(iUnit + unitCount, iCond) = pvals_per_cond(iCond);
                    end
                       %responsive cell?  shuffle pre/post labels 

                pre_activity = nanmean(curr_raster(:, 40:50), 2);
                post_activity = nanmean(curr_raster(:, 55:65), 2);
                all_activity = [pre_activity; post_activity];
                shuffling_idx = randi(size(all_activity,1),1000);
                shuffled_diffs = arrayfun(@(x) nanmean(all_activity(shuffling_idx(1:500,x))...
                    - all_activity(shuffling_idx(501:1000,x))), 1:1000);
                real_diff = nanmean(post_activity - pre_activity);
                pctile_2 = prctile(shuffled_diffs, 95);
                pctile_1 = prctile(shuffled_diffs, 5);

                passive_data_per_cond.pvalue_shuffled_005{experiment_type}(iUnit+unitCount) = ...
                   real_diff >= pctile_2 || real_diff <= pctile_1;

                 passive_data_per_cond.pvalue_shuffled_0025{experiment_type}(iUnit+unitCount) = ...
                   real_diff >= prctile(shuffled_diffs, 97.5) || real_diff <= prctile(shuffled_diffs, 2.5);

                                  colorMtx =  bc_colors(3, 'w');
%                 figure();
%                 subplot(6,2,[1,3, 5, 7])
%                 scatter(t(raster_x), raster_y, 2, [0,0,0], 'filled')
%                 ylabel('trial #')
%                 xticklabels({''})
%                 makepretty;
% 
%                 subplot(6,2,[9, 11])
%                 plot(t, curr_psth .* 1/psth_bin_size,'Color', colorMtx(2,:))
%                 xlabel('time from stim onset (s)')
%                 ylabel('spikes/s')
%                 legend(['sign rank p-value = ', num2str(passive_data_per_cond.pvalue{experiment_type}(iUnit+unitCount),3)])
%                 makepretty;
% 
%                 subplot(6,2,[6, 8])
%                 h = histogram(shuffled_diffs, 'FaceColor', colorMtx(1,:), 'EdgeColor', colorMtx(1,:), 'Normalization' ,'probability'); 
%                 hold on;
%                 ylims = ylim;
%                 line([real_diff, real_diff], [ylims(1), ylims(2)], 'Color', colorMtx(2,:))
%                 
%                 l1 = line([real_diff, real_diff], [ylims(1), ylims(2)], 'Color', colorMtx(2,:));
%                 l2 = line([pctile_2, pctile_2], [ylims(1), ylims(2)], 'Color', colorMtx(3,:));
%                 line([pctile_1, pctile_1], [ylims(1), ylims(2)], 'Color', colorMtx(3,:))
%                 xlabel('average f.r. post-stim - pre-stim')
%                 ylabel('fraction')
%                 legend([h, l1, l2], {'shuffled data', 'real value', '2.5%'})
%                 
% 
%                 makepretty;
% 

                    %psth

                    [curr_psth, ~, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(1:2:end), trial_cond_idx(1:2:end));

                    passive_data_per_cond.psth{experiment_type}(iUnit + unitCount, 1, :, :) = curr_psth;

                    [curr_psth, ~, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(2:2:end), trial_cond_idx(2:2:end));

                    passive_data_per_cond.psth{experiment_type}(iUnit + unitCount, 2, :, :) = curr_psth;
                    %disp(unique_templates(shank_units(units_to_keep(iUnit))))
                    
                    [curr_psth, curr_raster, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, ...
                        align_times(1:1:end), trial_cond_idx(1:1:end));
                     passive_data_per_cond.psth{experiment_type}(iUnit + unitCount, 3, :, :) = curr_psth;

                     psth_bin_size_det = 0.001;
                    raster_window_det = [-0.2, 0.6];
                    [curr_psth, curr_raster, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(shank_units(units_to_keep(iUnit))), raster_window_det, psth_bin_size_det, ...
                        align_times(1:1:end), trial_cond_idx(1:1:end));
                     startIdx = find(t >= 0.05, 1, 'first');
                     stopIdx = find(t >= 0.15, 1, 'first');
                     passive_data_per_cond.av_per_trial{experiment_type, iRecording}(iUnit,:) = nanmean(curr_raster(:,startIdx:stopIdx),2);
                     startIdx = find(t >= -0.15, 1, 'first');
                     stopIdx = find(t >= -0.05, 1, 'first');
                     passive_data_per_cond.av_per_trial_base{experiment_type, iRecording}(iUnit,:) = nanmean(curr_raster(:,startIdx:stopIdx),2);
                     passive_data_per_cond.av_psth{experiment_type, iRecording}(iUnit,:,:) = curr_psth;
                    
                     % area
                     %passive_data_per_cond.unit_area{experiment_type, iRecording}(iUnit+unitCount)
                     % unit type 

            
                end

            end

            % save data in structure
            if ~isempty(units_to_keep)
                passive_data_per_cond.nz_trial_types{experiment_type, iRecording} = trial_conditions;
                passive_data_per_cond.nz_no_move_trials{experiment_type, iRecording} = no_move_trials;
                passive_data_per_cond.animal_day_site_shank(unitCount+1:unitCount+size(units_to_keep, 1), :) = ...
                    repmat([mouse_day_sites_shank_rec(iRecording, 1), curr_day, site, curr_shank], size(units_to_keep, 1), 1);
                passive_data_per_cond.unit_area(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_area;
                passive_data_per_cond.unit_coords(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_coords;
                passive_data_per_cond.t = t;
                passive_data_per_cond.unitNum((unitCount + 1:unitCount + size(units_to_keep, 1))) = shank_units(units_to_keep);
                passive_data_per_cond.propISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                passive_data_per_cond.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression(units_to_keep);
                passive_data_per_cond.wvDur((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.templateDuration(units_to_keep);
                passive_data_per_cond.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.spike_rateSimple(units_to_keep);
                passive_data_per_cond.unitType((unitCount + 1:unitCount + size(units_to_keep, 1))) = unitType(units_to_keep);
                %passive_data_per_cond
                unitCount = unitCount + size(units_to_keep, 1);

            end

            % clear variables
            disp(['   ', num2str(iRecording), '/', num2str(length(use_recs))])

            keep mouse_day_sites_shank_rec unique_mice passive_data_per_cond use_recs info_table regions regions_id unitCount st load_type keep_type
        end
    end
    %catch
    %end
end

end