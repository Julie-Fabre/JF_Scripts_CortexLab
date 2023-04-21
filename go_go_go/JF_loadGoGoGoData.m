
%% load units in various regions + their quality metric clasification and thier location

myPaths;
[tv, av, st, bregma] = bd_loadAllenAtlas(atlasBrainRegLocation);

animals = {'JF096', 'JF097', 'JF099'};
regionsOfInterest = {'CP', 'GPe', 'SNr', 'GPi', 'STN'};
for iRegion = 1:length(regionsOfInterest)
    regionOfInterest_number(iRegion) = st.id(strcmp(st.acronym, regionsOfInterest{iRegion}));
end

protocols = {'noGo_path2stage3', 'JF_choiceworldStimuli_wheel'};

ephysData = struct;

for iAnimal = 1%:length(animals)
    keep animals iAnimal regionsOfInterest protocols ephysData regionOfInterest_number
    animal = animals{iAnimal};
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)

    experiments_active = AP_find_experimentsJF(animal, protocols{1}, true);
    experiments_active = experiments_active([experiments_active.ephys]);

    experiments_passive = AP_find_experimentsJF(animal, protocols{2}, true);
    experiments_passive = experiments_passive([experiments_passive.ephys]);
                
    iRec = 0;
    for iRegion = 1%:length(regionsOfInterest)
        
        for iProbe = 1:length(probe2ephys)

            this_ccf = probe_ccf(iProbe); %find if contains location of interest
            if any(ismember(this_ccf.trajectory_areas, regionOfInterest_number(iRegion)))
                iRec = iRec + 1;
                site = probe2ephys(iProbe).site;
                day = experiments_active(probe2ephys(iProbe).day).day;
                if isfield(probe2ephys, 'recording') && ~isempty(probe2ephys(iProbe).recording)
                    recording = probe2ephys(iProbe).recording(1);
                else
                    recording = [];
                end
                shank = probe2ephys(iProbe).shank;
                lfp = NaN;
                load_sync = false;
                loadClusters = 0;
                rerun = false;
                plotGUI = false;
                runQM = false;
    
                [unitType,qMetrics] = bc_qualityMetricsPipeline_JF(animal, day, site, 1, protocols{1}, rerun, plotGUI, runQM);
                
                %% active experiments
                % get experiments active dernier qui a au moins 100 trials
                experiment = experiments_active(probe2ephys(iProbe).day).experiment(end);
                JF_load_experiment;
                locationDepths = this_ccf.probe_depths(find(ismember(this_ccf.trajectory_areas, regionOfInterest_number(iRegion)))); %QQ change to detect if multiple chunks
                %theseDepths = [min(locationDepths), max(locationDepths)];
                theseUnits = template_depths >= min(locationDepths) & template_depths <= max(locationDepths);
                uniqueTemps = unique(spike_templates);
                theseSpikes = ismember(spike_templates, uniqueTemps(theseUnits));
                % if striatum, get prop_isi and waveform duration and
                % classify 
                theseUnits_idx = find(theseUnits);
                if iRegion == 1 % striatum 
                    for iUnit = 1:size(theseUnits_idx,1)
                        acg = CCGBz([spike_times_timeline(spike_templates==uniqueTemps(theseUnits_idx(iUnit)));...
                            spike_times_timeline(spike_templates==uniqueTemps(theseUnits_idx(iUnit)))],...
                            [ones(size(spike_times_timeline(spike_templates==uniqueTemps(theseUnits_idx(iUnit))),1),1);...
                            ones(size(spike_times_timeline(spike_templates==uniqueTemps(theseUnits_idx(iUnit))),1),1).*2],...
                            'binsize',0.001, 'duration', 1);
                        prop_isi(iUnit) = find(acg(500:end)>=nanmean(acg(600:900)),1,'first');
                        waveform_duration(iUnit) = qMetrics.waveformDuration_peakTrough(theseUnits_idx(iUnit));
                    end
                    figure()
                    scatter(prop_isi, waveform_duration)
                end

                ephysData(iRec).spike_times = spike_times_timeline(theseSpikes);
                ephysData(iRec).spike_templates = spike_templates(theseSpikes);
                ephysData(iRec).trial_conditions = trial_conditions;
                rewards_times = [signals_events.responseTimes(logical(abs(crrectExpect-1)))'; ... 
                    signals_events.responseTimes(crrectExpect)'; unexpected_reward_times'];% incorrect, correct, unexpected 
                rewards_type = [-1*ones(size(signals_events.responseTimes(logical(abs(crrectExpect-1))),2),1); ...
                    trial_outcome_reward_omission'; 2*ones(size(unexpected_reward_times,2),1)];
                ephysData(iRec).rewards_times = rewards_times;
                ephysData(iRec).rewards_type = rewards_type;
                ephysData(iRec).stimOn_times = stimOn_times;
                ephysData(iRec).stim_to_move = stim_to_move;
                ephysData(iRec).stim_to_feedback = stim_to_feedback;
                ephysData(iRec).unitType = unitType(theseUnits);
                ephysData(iRec).prop_isi =  prop_isi;
                ephysData(iRec).waveform_duration = waveform_duration;
                %% passive experiments
                % get experiments

                % loop through, load data, and save in structure
                experiment = experiments_passive(probe2ephys(iProbe).day).experiment(end);
                JF_load_experiment;
                locationDepths = this_ccf.probe_depths(find(ismember(this_ccf.trajectory_areas, regionOfInterest_number(iRegion)))); %QQ change to detect if multiple chunks
                %theseDepths = [min(locationDepths), max(locationDepths)];
                theseUnits = template_depths >= min(locationDepths) & template_depths <= max(locationDepths);
                uniqueTemps = unique(spike_templates);
                theseSpikes = ismember(spike_templates, uniqueTemps(theseUnits));

                ephysData(iRec).spike_times_passive = spike_times_timeline(theseSpikes);
                ephysData(iRec).spike_templates_passive = spike_templates(theseSpikes);
                ephysData(iRec).trial_conditions_passive = trial_conditions;
                ephysData(iRec).stimOn_times_passive = stimOn_times;
                ephysData(iRec).stim_to_move_passive = stim_to_move;
                

                ephysData(iRec).location = regionsOfInterest{iRegion};
            end
        end
    end
end