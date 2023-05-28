%% paramaters 
load_passive_ = true;
passive_info = readtable('/home/julie/Dropbox/PhD_summary/allPassiveRecs.csv');
regions = {'CP', 'GPe', 'GPi', 'STN', 'SNr', 'SNc', 'VTA'};
%regions = {'CP'};
warning off;
%get region nums (allen format)
cl_myPaths;
if ~exist('st', 'var')
[tv, av, st, bregma] = bd_loadAllenAtlas(atlasBrainRegLocation);
end
for iRegion = 1:size(regions,2)
    regions_id(iRegion) = st.id(strcmp(st.acronym, regions(iRegion)));
end

% get all mice x day x site combinations 
use_recs = strcmp(passive_info.Use_, 'Yes');

% good units (bombcell) (+ cell type?) . location. psth per image type.

%% load all mouse rec locations 
passive_info = sortrows(passive_info,1); % make sure data is sorted by mouse 
unique_mice = unique(passive_info.Mouse(use_recs));

% initialize 
mouse_day_sites_shank_rec = [];
disp('loading data...')
for iMouse = 1:length(unique_mice)
    day_sites_shank_rec = [passive_info.DayNum(use_recs & strcmp(passive_info.Mouse, unique_mice(iMouse))), ...
        passive_info.Site(use_recs & strcmp(passive_info.Mouse, unique_mice(iMouse))), ...
        passive_info.Shanks(use_recs & strcmp(passive_info.Mouse, unique_mice(iMouse))),...
        passive_info.Rec_Exp(use_recs & strcmp(passive_info.Mouse, unique_mice(iMouse)))];

    mouse_day_sites_shank_rec = [mouse_day_sites_shank_rec; repmat(iMouse, length(day_sites_shank_rec), 1),...
        day_sites_shank_rec];

end

%% load data 
passive_data_per_cond = struct;
passive_data_per_cond.psth = cell(5, 1);
passive_data_per_cond.pvalue = cell(5,1);
passive_data_per_cond.psth_conditions = cell(5, 1);
use_recs = find(use_recs);
unitCount = 0;
for iRecording = 1:length(use_recs)
    % curr variables 
    animal = unique_mice{mouse_day_sites_shank_rec(iRecording,1)};
    curr_day =  mouse_day_sites_shank_rec(iRecording,2);
    site = mouse_day_sites_shank_rec(iRecording,3);
    curr_shank =  mouse_day_sites_shank_rec(iRecording,4);
    curr_rec = mouse_day_sites_shank_rec(iRecording,5);

    % load probe_ccf and probe2ephys
    probe2ephys_location = AP_cortexlab_filenameJF(animal, [],[], 'probe2ephys');
    load(probe2ephys_location)
    day_sites_shank_probe = [probe2ephys.day; probe2ephys.site; probe2ephys.shank]';
    
    mouse_day_sites_shank_rec(isnan(mouse_day_sites_shank_rec(:,4)),4) =0 ;%replace nan by 0 
    day_sites_shank_probe(isnan(day_sites_shank_probe(:,3)),3) =0 ;%replace nan by 0 
    
    [probe_rec, probe_rec_idx] = ismember(mouse_day_sites_shank_rec(iRecording,2:4),day_sites_shank_probe, 'rows');

    probe_ccf_location = AP_cortexlab_filenameJF(animal, [],[], 'histo');
    load(probe_ccf_location)

    % is any region present 
    this_probe = probe_rec_idx;
    for iRegion = 1:size(regions,2)
        this_region = regions_id(iRegion);
        this_region_idx = find(probe_ccf(this_probe).trajectory_areas == this_region);
        if ~isempty(this_region_idx)
            this_region_start(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(1));
            this_region_stop(iRegion) = probe_ccf(this_probe).probe_depths(this_region_idx(end));
        else
            this_region_start(iRegion) = 0;
            this_region_stop(iRegion) = 0;
        end
    end
    
    % load units in relevant depths 
    these_regions_present = find(this_region_stop);
    if ~isempty(these_regions_present)
        these_exps = str2num(passive_info.Exps{use_recs(iRecording)});
        % in order: gratings, location, (exp type 1,2,3,4)
        % conditions = trial_conditions 
        for iExperiment = 1:size(these_exps,2)
            experiment = these_exps(iExperiment);%for now, just load 1 experiment
            experiments = AP_find_experimentsJF(animal, '', true);
            experiments = experiments([experiments.ephys]);
            day = experiments(curr_day).day;
            if isnan(curr_rec)
                recording = [];
            else
                recording = curr_rec;
            end
            
    
            verbose = false; % display load progress and some info figures
            load_parts.cam=false;
            load_parts.imaging=false;
            load_parts.ephys=true;
            JF_load_experiment;
            if contains(expDef, 'JF_GratingPassiveVarITI')
                experiment_type = 1;
            elseif contains(expDef, 'JF_locationsFit')
                experiment_type = 2;
            elseif contains(expDef, 'JF_natural_imagesFit')
                experiment_type = 3;
            elseif contains(expDef, 'JF_choiceworldStimuli_onlyTask')
                experiment_type = 5;
            elseif contains(expDef, 'JF_choiceworldStimuli')
                experiment_type = 4;
            end
            % "real" template depths 
    %         metaFile = strrep(ephysAP_path, '.cbin', '.meta');
    %         [~, channelMapIMRO] = bc_readSpikeGLXMetaFile(metaFile);
    %         if contains(channelMapIMRO, 'hStripe')
    %             template_depths = template_depths - 2880 + 750;
    %         end
    
            %subselect shank 
            if iExperiment == 1
                if curr_shank > 0
                    shank_xdepth = [250 * (curr_shank-1), 250 * (curr_shank-1) + 32];
                    shank_units = find(template_xdepths >= shank_xdepth(1) & template_xdepths <= shank_xdepth(2));
                else
                    shank_units = 1:size(template_xdepths,1)';
                end
                % get units we want to keep + store each units location (region + coordinates)
                units_to_keep = [];
                units_to_keep_area = [];
                units_to_keep_coords = [];
                for iRegion = 1:length(these_regions_present)
                    new_units =  find(template_depths(shank_units) >= this_region_start(these_regions_present(iRegion)) &...
                        template_depths(shank_units) <= this_region_stop(these_regions_present(iRegion)) );
                    units_to_keep = [units_to_keep; new_units];
                    units_to_keep_area = [units_to_keep_area; ones(size(new_units,1),1).*these_regions_present(iRegion)];
                    
                    unit_closest_depth = arrayfun(@(x) ...%closest depth 
                       find(probe_ccf(this_probe).probe_depths >= template_depths(shank_units(new_units(x))), 1, 'first'), 1:length(new_units));
                    units_to_keep_coords = [units_to_keep_coords; ...
                        probe_ccf(this_probe).trajectory_coords(unit_closest_depth,:)];
                end
            end

        
            % get stim response 
            unique_templates = unique(spike_templates);
            passive_data_per_cond.psth_conditions{experiment_type} = unique(trial_conditions, 'rows');
            % for half trials, identify max. 
            % then average on that half. 
            for iUnit = 1:size(units_to_keep,1)
                raster_window = [-0.5, 1];
                psth_bin_size = 0.01;
                align_times = stimOn_times;
                [~, trial_cond_idx] = ismember(trial_conditions,  passive_data_per_cond.psth_conditions{experiment_type}, 'rows');
                %responsive cell 

                
                 [~, curr_raster, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline,...
                    unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size,...
                    align_times, []);
                 passive_data_per_cond.pvalue{experiment_type}(iUnit + unitCount) =...
                     signrank(nanmean(curr_raster(:,20:40),2), nanmean(curr_raster(:,55:75),2));
                 


                %psth
                
                [curr_psth, ~, ~, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline,...
                    unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size,...
                    align_times(1:2:end), trial_cond_idx(1:2:end));

                passive_data_per_cond.psth{experiment_type}(iUnit + unitCount, 1, :, :) = curr_psth;

                [curr_psth, ~, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline,...
                    unique_templates(shank_units(units_to_keep(iUnit))), raster_window, psth_bin_size,...
                    align_times(2:2:end), trial_cond_idx(2:2:end));

                passive_data_per_cond.psth{experiment_type}(iUnit + unitCount, 2, :, :) = curr_psth;

                


            end
        end

    % save data in structure    
    if ~isempty(new_units)
        passive_data_per_cond.animal_day_site_shank(unitCount+1:unitCount+size(units_to_keep,1),:) = ...
            repmat([mouse_day_sites_shank_rec(iRecording,1), curr_day, site, curr_shank], size(units_to_keep,1), 1);
        passive_data_per_cond.unit_area(unitCount+1:unitCount+size(units_to_keep,1),:) = units_to_keep_area;
        passive_data_per_cond.unit_coords(unitCount+1:unitCount+size(units_to_keep,1),:) = units_to_keep_coords;
        passive_data_per_cond.t = t;
        passive_data_per_cond.unitNum((unitCount+1:unitCount+size(units_to_keep,1))) = shank_units(units_to_keep);
        unitCount = unitCount + size(units_to_keep,1);
    
    end

    % clear variables 
    disp(['   ', num2str(iRecording), '/', num2str(length(use_recs))])
    
    keep mouse_day_sites_shank_rec unique_mice passive_data_per_cond use_recs passive_info regions regions_id unitCount st
    end
end