function passive_data = cl_loadAverageData(load_type)

%% paramaters
load_passive = true;
runEP = true;
if contains(load_type, 'passive')
    info_table = readtable('/home/julie/Dropbox/PhD_summary/allPassiveRecs.csv', 'VariableNamingRule', 'modify');
elseif contains(load_type, 'task')
    info_table = readtable('/home/julie/Dropbox/PhD_summary/allTaskRecs.csv', 'VariableNamingRule', 'modify');
end
regions = {'CP', 'GPe', 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
%regions = {'CP'};

%get region nums (allen format)
cl_myPaths;
if ~exist('st', 'var')
    [~, ~, st, ~] = bd_loadAllenAtlas(atlasBrainRegLocation);
end
for iRegion = 1:size(regions, 2)
    regions_id(iRegion) = st.id(strcmp(st.acronym, regions(iRegion)));
end

% get all mice x day x site combinations
use_recs = strcmp(info_table.Use_, 'Yes');

% good units (bombcell) (+ cell type?) . location. psth per image type.

%% load all mouse rec locations
info_table = sortrows(info_table, 1); % make sure data is sorted by mouse
unique_mice = unique(info_table.Mouse(use_recs));

% initialize
mouse_day_sites_shank_rec = [];

for iMouse = 1:length(unique_mice)
    day_sites_shank_rec = [info_table.DayNum(use_recs & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Site(use_recs & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Shanks(use_recs & strcmp(info_table.Mouse, unique_mice(iMouse))), ...
        info_table.Rec_Exp(use_recs & strcmp(info_table.Mouse, unique_mice(iMouse)))];

    mouse_day_sites_shank_rec = [mouse_day_sites_shank_rec; repmat(iMouse, size(day_sites_shank_rec, 1), 1), ...
        day_sites_shank_rec];

end

%% load data
passive_data = struct;
use_recs = find(use_recs);
unitCount = 0;

for iRecording = 1:length(use_recs)
    % curr variables
    animal = unique_mice{mouse_day_sites_shank_rec(iRecording, 1)};
    curr_day = mouse_day_sites_shank_rec(iRecording, 2);
    site = mouse_day_sites_shank_rec(iRecording, 3);
    curr_shank = mouse_day_sites_shank_rec(iRecording, 4);
    curr_rec = mouse_day_sites_shank_rec(iRecording, 5);

    % load probe_ccf and probe2ephys
    probe2ephys_location = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys');
    load(probe2ephys_location)
    day_sites_shank_probe = [probe2ephys.day; probe2ephys.site; probe2ephys.shank]';

    mouse_day_sites_shank_rec(isnan(mouse_day_sites_shank_rec(:, 4)), 4) = 0; %replace nan by 0
    day_sites_shank_probe(isnan(day_sites_shank_probe(:, 3)), 3) = 0; %replace nan by 0

    [probe_rec, probe_rec_idx] = ismember(mouse_day_sites_shank_rec(iRecording, 2:4), day_sites_shank_probe, 'rows');

    probe_ccf_location = AP_cortexlab_filenameJF(animal, [], [], 'histo');
    load(probe_ccf_location)

    % is any region present
    this_probe = probe_rec_idx;
    if this_probe >= 1
        if isfield(probe_ccf, 'probe_depths')
            for iRegion = 1:size(regions, 2)
                this_region = regions_id(iRegion);
                if ~isempty(probe_ccf(this_probe).trajectory_areas)
                    this_region_idx = find(probe_ccf(this_probe).trajectory_areas == this_region);
                    if ~isempty(this_region_idx) && ~isempty(probe_ccf(this_probe).probe_depths)
                        this_region_start(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(1));
                        this_region_stop(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(end));
                    else
                        if isfield('probe_depths', probe_ccf)
                            if isempty(probe_ccf(this_probe).probe_depths)
                                warning on;
                                warning(['probe_ccf depths empty: ', animal, ', day: ', ...
                                    num2str(curr_day), ', site: ', num2str(site), ', probe: ', num2str(this_probe)])
                                warning off;
                            end
                        end
                        this_region_start(iRegion) = 0;
                        this_region_stop(iRegion) = 0;
                    end
                end
            end
        else

            this_region_start(1:size(regions, 2)) = 0;
            this_region_stop(1:size(regions, 2)) = 0;
        end

        % load units in relevant depths
        these_regions_present = find(this_region_stop);
        if ~isempty(these_regions_present)
            these_exps = str2num(info_table.Exps{use_recs(iRecording)});
            experiment = these_exps(1); %for now, just load 1 experiment
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
            loadClusters = false; %QQ for now. add bombcell output later
            JF_load_experiment;
            % "real" template depths
            %         metaFile = strrep(ephysAP_path, '.cbin', '.meta');
            %         [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);
            %         if contains(channelMapIMRO, 'hStripe')
            %             template_depths = template_depths - 2880 + 750;
            %         end

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


            % get stim response
            unique_templates = unique(spike_templates);

            for iUnit = 1:size(units_to_keep, 1)
                raster_window = [-0.5, 1];
                psth_bin_size = 0.01;
                align_times = stimOn_times;
                [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                    unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size, align_times, []);
                
                passive_data.psth(iUnit+unitCount, :) = curr_psth;

                %responsive cell? signrank test
                 passive_data.pvalue(iUnit+unitCount) = ...
                    signrank(nanmean(curr_raster(:, 40:50), 2), nanmean(curr_raster(:, 55:75), 2));
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

                passive_data.pvalue_shuffled_005(iUnit+unitCount) = ...
                   real_diff >= pctile_2 || real_diff <= pctile_1;

                 passive_data.pvalue_shuffled_0025(iUnit+unitCount) = ...
                   real_diff >= prctile(shuffled_diffs, 97.5) || real_diff <= prctile(shuffled_diffs, 2.5);
                
%                 colorMtx =  bc_colors(3, 'w');
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
%                 legend(['sign rank p-value = ', num2str(passive_data.pvalue(iUnit+unitCount),3)])
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
%                 %shuffling_idx = arrayfun(@(x) randperm(size(pre_activity,1)), 1:1000);

               
            end
            % get ephys prop: firing rate, prop ISI, pss, template dur
            protocol = '';
            rerunEP = 0;
            plotGUI = 0;
            runQM = 1;
            rerunQM = 0;
            region = '';
            try
                [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, day, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
            catch
                rerunQM = 1;
                [unitType, qMetrics] = bc_qualityMetricsPipeline_JF(animal, day, site, recording, 1, protocol, rerunQM, plotGUI, runQM);
            end
            %         warning('no quality metrics ')
            %         unitType = nan(size(unique_templates,1), 1);
            %     end

            ephysProperties = bc_ephysPropertiesPipeline_JF(animal, day, site, recording, 1, protocol, rerunEP, runEP, region);
            % save data in structure
            if ~isempty(new_units)
                passive_data.animal_day_site_shank(unitCount+1:unitCount+size(units_to_keep, 1), :) = ...
                    repmat([mouse_day_sites_shank_rec(iRecording, 1), curr_day, site, curr_shank], size(units_to_keep, 1), 1);
                passive_data.unit_area(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_area;
                passive_data.unit_coords(unitCount+1:unitCount+size(units_to_keep, 1), :) = units_to_keep_coords;
                passive_data.t = t;
                passive_data.unitNum((unitCount + 1:unitCount + size(units_to_keep, 1))) = shank_units(units_to_keep);
                passive_data.propISI((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.propLongISI(units_to_keep);
                passive_data.pss((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.postSpikeSuppression(units_to_keep);
                passive_data.wvDur((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.templateDuration(units_to_keep);
                passive_data.fr((unitCount + 1:unitCount + size(units_to_keep, 1))) = ephysProperties.spike_rateSimple(units_to_keep);
                passive_data.unitType((unitCount + 1:unitCount + size(units_to_keep, 1))) = unitType(units_to_keep);
                unitCount = unitCount + size(units_to_keep, 1);

            end
        end

        % clear variables
        keep mouse_day_sites_shank_rec unique_mice passive_data use_recs info_table regions regions_id unitCount st runEP iRecording
    end
    disp(num2str(iRecording))
end

end