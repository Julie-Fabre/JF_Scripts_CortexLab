function allDataStruct = JF_loadAllData(animals, region, protocol, useQualityMetrics)
myPaths;
allDataStruct  = struct;
% TO DO:
% dealing with 4 shanks
% using old brains st things
% channel locations 
for iAnimal = 1:size(animals)

    animal = animals{iAnimal};

    % find experiments for this animal
    experiments = AP_find_experimentsJF(animal, protocol, true);
    experiments = experiments([experiments.ephys]);

    % get histology information
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)

    [~, ~, st, ~] = JF_loadAllenAtlasFiles(probe2ephysFile);

    for iRecording = 1:size(probe2ephys, 2)
        % get recording info
        curr_day = probe2ephys(iRecording).day; % (set which day to use)
        day = experiments(curr_day).day; % (set which day to use)
        
        site = probe2ephys(iRecording).site; % (set which site to use)
        curr_shank = probe2ephys(iRecording).shank; % (set which day to use)
        if isfield(probe2ephys, 'recording')
            recording = probe2ephys(iRecording).recording; % (set which day to use)
        else
            recording = [];
        end

        % check histology is aligned
        if isempty(probe_ccf(iRecording).probe_depths)
            disp('histology not aligned, skipping this experiment')
            continue;
        end

        % check region of interest is in this experiment

        ephysDir = AP_cortexlab_filenameJF(animal, day, [], 'ephys', site, recording);
        channel_positions = readNPY([ephysDir, filesep, 'channel_positions.npy']);
        topLimit = find(probe_ccf(iRecording).probe_depths >= min(channel_positions(:, 2)), 1, 'first');
        bottomLimit = find(probe_ccf(iRecording).probe_depths >= max(channel_positions(:, 2)), 1, 'first');
        region_id = st.id(find(~cellfun(@isempty, strfind(st.acronym, region))));
        area_contained = find(probe_ccf(iRecording).trajectory_areas == region_id);
        locationKeep = probe_ccf(iRecording).probe_depths(area_contained);

        if ~any(area_contained >= topLimit & area_contained <= bottomLimit)
            disp('region of interest not present, skipping this experiment')
            continue;
        end

        % loading parameters
        verbose = false; % display load progress and some info figures
        load_parts.cam = false;
        load_parts.imaging = false;
        load_parts.ephys = true;

        % which experiment # to use
        n_trials = zeros(size(experiments(curr_day).experiment, 2), 1);
        for iExperiment = experiments(curr_day).experiment
            [block_filename, ~] = AP_cortexlab_filenameJF(animal, day, iExperiment, 'block');
            try
                load(block_filename)
                if isfield(block.events, 'stim_idValues')
                    n_trials(iExperiment) = length(block.events.stim_idValues);
                elseif isfield(block.events, 'stimulusOnTimes')
                    n_trials(iExperiment) = length(block.events.stimulusOnTimes);
                end
            catch
                n_trials(iExperiment) = NaN;
            end
        end
        experiment = find(n_trials == max(n_trials));

        % load quality metrics if present
        ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
        savePath = fullfile(ephysDirPath, 'qMetrics');
        qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
        if ~isempty(qMetricsExist) && useQualityMetrics
            qMetric = parquetread(fullfile(savePath, 'templates._bc_qMetrics.parquet'));
            param = parquetread(fullfile(savePath, '_bc_parameters._bc_qMetrics.parquet'));
            bc_getQualityUnitType;
        end

        % load experiment
        JF_load_experiment;

        % if shank, keep only shank data 
        if ~isnan(curr_shank)
            theseChannelPositions = [(curr_shank-1) * 250, (curr_shank-1)*250 + 32];
            theseChannels = ismember(channel_positions(:,1), theseChannelPositions);
            theseTemplates = ismember(template_xdepths, theseChannelPositions);
            theseSpikes = ismember(spike_xdepths, theseChannelPositions);
            spike_times = spike_times(theseSpikes);
            spike_templates = spike_templates(theseSpikes);
            %rename 
            good_templates_idx = unique(spike_templates);
            new_spike_idx = nan(max(spike_templates), 1);
            new_spike_idx(good_templates_idx) = 1:length(good_templates_idx);
            spike_templates = new_spike_idx(spike_templates);
            
            template_amplitudes = template_amplitudes(theseSpikes);
            template_depths = template_depths(theseTemplates);
            channel_positions = channel_positions(theseChannels,:);
            templates = templates(theseTemplates,:,theseChannels); 
        end

        % get each unit s location 


        % store loaded data
        allDataStruct(iRecording).spike_templates = spike_templates;
        allDataStruct(iRecording).spike_times = spike_times_timeline;
        %allDataStruct(iRecording).spike_waveforms = templates
        %allDataStruct(iRecording).spike_locations
        %allDataStruct(iRecording).stimIDs 
        allDataStruct(iRecording).protocol = expDef;
        allDataStruct(iRecording).recording = iRecording;
        allDataStruct(iRecording).trial_conditions = trial_conditions;
        allDataStruct(iRecording).stimOn_times = stimOn_times;
        allDataStruct(iRecording).stim_to_move = stim_to_move;


    end

end

end