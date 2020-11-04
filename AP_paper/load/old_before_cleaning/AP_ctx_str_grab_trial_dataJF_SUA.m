
%% ~~~ AP_ctx_str_grab_trial_data: pull out data organized by trials ~~~

%% Set flags
trial_data = struct;
% filter_mua: convolve MUA with filter to match widefield fluorescence
if ~exist('filter_mua', 'var')
    filter_mua = false;
end

%% Get if task if passive dataset

% Task dataset if signals (expDef) and contains 'vanillaChoiceworld'
task_dataset = exist('expDef', 'var') && contains(expDef, 'vanillaChoiceworld');

%% Set parameters for cortical fluoresence

if verbose
    disp('Re-casting and deconvolving fluorescence...');
end

% Convert U to master U
load('E:\for_Julie\wf_processing\wf_alignment\U_master.mat');
Udf_aligned = AP_align_widefieldJF(Udf, animal, day);
fVdf_recast = ChangeU(Udf_aligned, fVdf, U_master);

% Set components to keep
use_components = 1:200;

% Deconvolve fluoresence (for use in regression)
fVdf_deconv = AP_deconv_wfJF(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;

%% Set parameters for striatal multiunit

n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

%% Set parameters for trials

% Get event-aligned activity
raster_window = [-0.5, 2];
upsample_factor = 1;
raster_sample_rate = 1 / (framerate * upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;


t_peri_event = bsxfun(@plus, use_align, t);
t_peri_event_bins = [t_peri_event - raster_sample_rate / 2, t_peri_event(:, end) + raster_sample_rate / 2];

% Pick trials to keep
if task_dataset
    use_trials = ...
        trial_outcome ~= 0 & ...
        ~signals_events.repeatTrialValues(1:n_trials)' & ...
        stim_to_feedback < 1.5;
    use_align_move = wheel_move_time;
use_align_move(isnan(use_align_move)) = 0;
use_align_outcome = signals_events.responseTimes';
use_align_outcome(isnan(use_align_outcome)) = 0;
t_peri_event_move = bsxfun(@plus, use_align_move, t);
t_peri_event_bins_move = [t_peri_event_move - raster_sample_rate / 2, t_peri_event_move(:, end) + raster_sample_rate / 2];
t_peri_event_outcome = bsxfun(@plus, use_align_outcome, t);
t_peri_event_bins_outcome = [t_peri_event_outcome - raster_sample_rate / 2, t_peri_event_outcome(:, end) + raster_sample_rate / 2];

else
    use_trials = true(size(stimIDs));
end

%% Get trial info

trial_info = struct;

if task_dataset

    stim_contrastside = signals_events.trialSideValues(1:n_trials)' .* ...
        signals_events.trialContrastValues(1:n_trials)';
    trial_info.stimulus = stim_contrastside(use_trials);

    trial_info.response = 3 - (abs((trial_choice(use_trials) + 1)/2) + 1);
    trial_info.repeatNum = ones(sum(use_trials), 1);

    trial_info.outcome = reshape(trial_outcome(use_trials), [], 1);

    trial_info.stim_to_move = stim_to_move(use_trials);

else % If passive dataset

    trial_info.stimulus = stimIDs;

end

%% Load cell-types
previous = 0;
if ~previous
    celltype_path = fullfile('E:/', animal, thisDate, '/analysis/');

    if exist([celltype_path, 'qMetric.mat'], 'file') && redo == 0
        load([celltype_path, 'qMetric.mat']);
        load([celltype_path, 'ephysParams.mat']);
        load([celltype_path, 'param.mat']);
        load([celltype_path, 'ephysData.mat']);
    else
        thisAnimal = animal;
        thisExperiment = experiment;
        exampleQualityMetCellType_single;
    end

    goodUnits = qMetric.numSpikes >= param.minNumSpikes & qMetric.waveformRawAmpli .* 0.195 >= param.minAmpli & ...
        qMetric.fractionRPVchunk <= param.maxRPV  ...
        & ephysParams.somatic == param.somaCluster; %waveforms, amplis, isolation, isiviolations
    


    %QQ add isoD metrics 
  
    % goodUnits = qMetric.numSpikes >= param.minNumSpikes & qMetric.waveformRawAmpli .* 0.195 >= param.minAmpli & ...
    %     qMetric.spatialDecayTemp1 >= param.minSpatDeKlowbound & qMetric.fractionRPVchunk <= param.maxRPV & ...
    %     qMetric.numPeaksTroughsTemp <= param.maxNumPeak & ...
    %     ephysParams.somatic == param.somaCluster & ephysParams.templateDuration < 800;
        
    msn = ephysParams.postSpikeSuppressionBf < 40 & ephysParams.templateDuration > param.cellTypeDuration;
    fsi = ephysParams.prop_long_isi < 0.15 & ...
        ephysParams.postSpikeSuppressionBf < param.cellTypePostS & ...
        ephysParams.templateDuration < param.cellTypeDuration;
    tan = ephysParams.postSpikeSuppressionBf > param.cellTypePostS & ephysParams.templateDuration > param.cellTypeDuration;
    tan2 = ephysParams.postSpikeSuppressionBf >= 40 & ephysParams.postSpikeSuppressionBf <= param.cellTypePostS & ephysParams.templateDuration > param.cellTypeDuration ;
    th = ephysParams.prop_long_isi < 0.8 & ephysParams.prop_long_isi >0.15 &...
        ephysParams.postSpikeSuppressionBf < param.cellTypePostS & ...
        ephysParams.templateDuration < param.cellTypeDuration;
    noise = ephysParams.prop_long_isi > 0.8 &...
        ephysParams.postSpikeSuppressionBf < param.cellTypePostS & ...
        ephysParams.templateDuration < param.cellTypeDuration;
    
 
    n_celltypes = 6;

    if param.strOnly
        theseMSNs_strIdx = find(msn);
        theseFSIs_strIdx = find(fsi);
        theseTANs_strIdx = find(tan);
        theseMSNs = qMetric.thisUnit(theseMSNs_strIdx);
        theseFSIs = qMetric.thisUnit(theseFSIs_strIdx);
        theseTAN2s = qMetric.thisUnit(theseFSIs_strIdx);
    else
        theseMSNs = find(msn);
        theseFSIs = find(fsi);
        theseTANs = find(tan);
        theseTAN2s = find(tan2);
        theseTHs = find(th);
        theseNoises = find(noise);
    end

    % get all ctype spike idxs

    %%QQ add counts
    theseNeurons = [theseMSNs, theseFSIs, theseTANs, theseTHs, theseNoises, theseTAN2s];
    
 
    %checking 
    %trial_data.postSpikeSuppressionBf(size(theseMSNs,2)+size(theseFSIs,2)+1:size(theseMSNs,2)+size(theseFSIs,2)+size(theseTANs,2))
    
    
    allGroupsSpikes = nan(size(spike_templates, 1), 1);
    allGroups = nan(numel(theseNeurons), 1);
    allAnimals = nan(numel(theseNeurons), 1);
    allDepths = nan(numel(theseNeurons), 1);
    allDays = nan(numel(theseNeurons), 1);
    allExperiments = nan(numel(theseNeurons), 1);
    %%QQ add counts
    msnCount = 0;
    theseNeurons = [];
    for iMSN = 1:size(theseMSNs, 2)
        theseSpikes = spike_templates == theseMSNs(iMSN);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)

            msnCount = msnCount + 1;
            allGroups(msnCount) = thisDepth;
            allDepths(msnCount) = template_depths(theseMSNs(iMSN));
            allAnimals(msnCount)=curr_animal;
            allDays(msnCount)=curr_day;
            allExperiments(msnCount)=experiment;
            theseNeurons = [theseNeurons, theseMSNs(iMSN)];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth; %1 - 4:MSNs
            end

        end
    end
    fsiCount = 0;
    for iFSI = 1:size(theseFSIs, 2)
        theseSpikes = spike_templates == theseFSIs(iFSI);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            fsiCount = fsiCount + 1;
            allGroups(msnCount+fsiCount) = thisDepth + n_depths;
            allDepths(msnCount+fsiCount) = template_depths(theseFSIs(iFSI));
            allAnimals(msnCount+fsiCount)=curr_animal;
            allDays(msnCount+fsiCount)=curr_day;
            allExperiments(msnCount+fsiCount)=experiment;
            theseNeurons = [theseNeurons, theseFSIs(iFSI)];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + n_depths; %5 - 8:FSIs
            end

        end
    end

    tanCount = 0;
    for iTAN = 1:size(theseTANs, 2)
        theseSpikes = spike_templates == theseTANs(iTAN);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            tanCount = tanCount + 1;
            allGroups(msnCount+fsiCount+tanCount) = thisDepth + (n_depths * 2);
            allDepths(msnCount+fsiCount+tanCount) = template_depths(theseTANs(iTAN));
            allAnimals(msnCount+fsiCount+tanCount)=curr_animal;
            allDays(msnCount+fsiCount+tanCount)=curr_day;
            allExperiments(msnCount+fsiCount+tanCount)=experiment;
            theseNeurons = [theseNeurons, theseTANs(iTAN)];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + (n_depths * 2); %9 - 12:TANs
            end

        end
    end
    
    thCount = 0;
    for iTH = 1:size(theseTHs, 2)
        theseSpikes = spike_templates == theseTHs(iTH);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            thCount = thCount + 1;
            allGroups(msnCount+fsiCount+tanCount+thCount) = thisDepth + (n_depths * 3);
            allDepths(msnCount+fsiCount+tanCount+thCount) = template_depths(theseTHs(iTH));
            allAnimals(msnCount+fsiCount+tanCount+thCount)=curr_animal;
            allDays(msnCount+fsiCount+tanCount+thCount)=curr_day;
            allExperiments(msnCount+fsiCount+tanCount+thCount)=experiment;
            theseNeurons = [theseNeurons, theseTHs(iTH)];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + (n_depths * 3); %9 - 12:TANs
            end

        end
    end
    
  
    noiseCount = 0;
    for iNoise = 1:size(theseNoises, 2)
        theseSpikes = spike_templates == theseNoises(iNoise);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            noiseCount = noiseCount + 1;
            allGroups(msnCount+fsiCount+tanCount +thCount +noiseCount) = thisDepth + ...
                (n_depths * 4);
            allDepths(msnCount+fsiCount+tanCount+thCount+noiseCount) =...
                template_depths(theseNoises(iNoise));
            allAnimals(msnCount+fsiCount+tanCount+thCount+noiseCount)=curr_animal;
            allDays(msnCount+fsiCount+tanCount+thCount+noiseCount)=curr_day;
            allExperiments(msnCount+fsiCount+tanCount+thCount+noiseCount)=experiment;
            theseNeurons = [theseNeurons, theseNoises(iNoise)];
            
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + (n_depths * 4); %9 - 12:TANs
            end

        end
    end
     
    tan2Count = 0;
    for iTAN2 = 1:size(theseTAN2s, 2)
        theseSpikes = spike_templates == theseTAN2s(iTAN2);

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            tan2Count = tan2Count + 1;
            allGroups(msnCount+fsiCount+tanCount+thCount+noiseCount+tan2Count) = thisDepth + (n_depths * 5);
            allDepths(msnCount+fsiCount+tanCount+thCount+noiseCount+tan2Count) = template_depths(theseTAN2s(iTAN2));
            allAnimals(msnCount+fsiCount+tanCount+thCount+noiseCount+tan2Count)=curr_animal;
            allDays(msnCount+fsiCount+tanCount+thCount+noiseCount+tan2Count)=curr_day;
            allExperiments(msnCount+fsiCount+tanCount+thCount+noiseCount+tan2Count)=experiment;
            theseNeurons = [theseNeurons, theseTAN2s(iTAN2)];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + (n_depths * 5); %9 - 12:TANs
            end

        end
    end     
    
    trial_data.numSpikes = qMetric.numSpikes(theseNeurons);
    trial_data.waveformRawAmpli = qMetric.waveformRawAmpli(theseNeurons);
    trial_data.fractionRPVchunk = qMetric.fractionRPVchunk(theseNeurons);
    trial_data.somatic = ephysParams.somatic(theseNeurons);
    trial_data.spatialDecayTemp1 = qMetric.spatialDecayTemp1(theseNeurons);
    trial_data.templateDuration = ephysParams.templateDuration(theseNeurons);
    trial_data.useTimeChunk = qMetric.useTimeChunk(theseNeurons,:);
    trial_data.symSpikesMissing = qMetric.symSpikesMissing(theseNeurons);
    trial_data.percent_missing_ndtr = qMetric.percent_missing_ndtr(:,theseNeurons);
    trial_data.timeChunks = ephysData.timeChunks;
    trial_data.numPeaksTroughs = qMetric.numPeaksTroughs(theseNeurons);
    trial_data.acg = ephysParams.ACG(theseNeurons,:); 
    trial_data.wv = qMetric.waveformUnit(theseNeurons,:);
    trial_data.goodUnits = goodUnits(theseNeurons); 
    trial_data.postSpikeSuppressionBf = ephysParams.postSpikeSuppressionBf(theseNeurons);
    trial_data.prop_long_isi = ephysParams.prop_long_isi(theseNeurons);
    
    
%     ff=find(allGroups== 9);
%     figure();plot(nanmean(trial_data.acg(ff,:)))
%     figure();plot(nanmean(ephysParams.ACG(theseTANs,:)))
else

    goodUnits = qMetrics(thisCount).numSpike > 300 & qMetricsAddI(thisCount).uQ > 1 & qMetricsAdd(thisCount).nmP < 4 &...
        qMetrics(thisCount).waveformRawAmpli > 80 ...
        & qMetrics(thisCount).fractionRPVchunk < 0.05 ...
        & qMetrics(thisCount).somaCluster == 1 & min(qMetrics(thisCount).percMssgAmpli(:, :)') < 30; %waveforms, amplis, isolation, isiviolations

     
    AlldurationC = qMetrics(thisCount).duration;
    AllpSC = ephysParams(thisCount).postSpikeSuppression;
    durLim = 400;
    maxLim=50;
    theseFSIs = find(goodUnits & AlldurationC < durLim) ;
    theseTANs = find(goodUnits & AlldurationC > durLim & AllpSC > maxLim);
    theseMSNs =  find(goodUnits & AlldurationC > durLim & AllpSC < maxLim);
    % get all ctype spike idxs

    %%QQ add counts
    theseNeurons = [theseMSNs, theseFSIs, theseTANs];

    allGroupsSpikes = nan(size(spike_templates, 1), 1);
    allGroups = nan(numel(theseNeurons), 1);
    allAnimals = nan(numel(theseNeurons), 1);
    allDepths = nan(numel(theseNeurons), 1);
    allDays = nan(numel(theseNeurons), 1);
    allExperiments = nan(numel(theseNeurons), 1);
    %%QQ add counts
    msnCount = 0;
    theseNeurons = [];
    [strUnits, ~] = find(ephysData(thisCount).str_templates);
    for iMSN = 1:size(theseMSNs, 2)
        
        %% load data 
        thisUnit = strUnits(theseMSNs(iMSN));
        theseSpikes = spike_templates == thisUnit;

        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)

            msnCount = msnCount + 1;
            allGroups(msnCount) = thisDepth;
            allDepths(msnCount+fsiCount) = template_depths(theseFSIs(iFSI));
            allAnimals(msnCount+fsiCount)=animal;
            allDays(msnCount+fsiCount)=thisDate;
            allExperiments(msnCount+fsiCount)=thisDate;
            theseNeurons = [theseNeurons, thisUnit];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth; %1 - 4:MSNs
                
            end

        end
    end
    fsiCount = 0;
    for iFSI = 1:size(theseFSIs, 2)
         thisUnit = strUnits(theseFSIs(iFSI));
        theseSpikes = spike_templates == thisUnit;


        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            fsiCount = fsiCount + 1;
            allGroups(msnCount+fsiCount) = thisDepth + n_depths;
             allDepths(msnCount+fsiCount) = template_depths(theseFSIs(iFSI));
            allAnimals(msnCount+fsiCount)=animal;
            allDays(msnCount+fsiCount)=thisDate;
            allExperiments(msnCount+fsiCount)=thisDate;
            theseNeurons = [theseNeurons, thisUnit];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + n_depths; %5 - 8:FSIs
            end

        end
    end

    tanCount = 0;
    for iTAN = 1:size(theseTANs, 2)
         thisUnit = strUnits(theseTANs(iTAN));
        theseSpikes = spike_templates == thisUnit;


        %get depth group
        thisDepth = unique(depth_group(theseSpikes));
        if ~isnan(thisDepth)
            tanCount = tanCount + 1;
            allGroups(msnCount+fsiCount+tanCount) = thisDepth + (n_depths * 2);
            allDepths(msnCount+fsiCount+tanCount) = template_depths(theseTANs(iTAN));
            allAnimals(msnCount+fsiCount+tanCount)=animal;
            allDays(msnCount+fsiCount+tanCount)=thisDate;
            allExperiments(msnCount+fsiCount+tanCount)=thisDate;
            theseNeurons = [theseNeurons, thisUnit];
            if numel(thisDepth) > 1
                error('more than one depth group for this unit. something is wrong')
            else
                allGroupsSpikes(theseSpikes) = thisDepth + (n_depths * 2); %9 - 12:TANs
            end

        end
    end
end

%% Trial-align data

if verbose;
    disp('Trial-aligning data...');
end;

% Cortical fluorescence
event_aligned_V = ...
    interp1(frame_t, fVdf_recast(use_components, :)', t_peri_event);

% Striatal multiunit
event_aligned_mua = nan(length(stimOn_times), length(t), numel(theseNeurons));

for iDepth = 1:numel(theseNeurons)
    thisNeuron = theseNeurons(iDepth);
    curr_spikes = spike_times_timeline(spike_templates == thisNeuron);

    % (skip if no spikes at this depth)
    if isempty(curr_spikes)
        %disp('empty')
        continue
    end

    event_aligned_mua(:, :, iDepth) = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes, t_peri_event_bins(x, :)), ...
        [1:size(t_peri_event, 1)]', 'uni', false)) ./ raster_sample_rate;
%     event_aligned_move_mua(:, :, iDepth) = cell2mat(arrayfun(@(x) ...
%         histcounts(curr_spikes, t_peri_event_bins_move(x, :)), ...
%         [1:size(t_peri_event_move, 1)]', 'uni', false)) ./ raster_sample_rate;
%     event_aligned_outcome_mua(:, :, iDepth) = cell2mat(arrayfun(@(x) ...
%         histcounts(curr_spikes, t_peri_event_bins_outcome(x, :)), ...
%         [1:size(t_peri_event_outcome, 1)]', 'uni', false)) ./ raster_sample_rate;

end

%
% (filter MUA if selected)
if filter_mua
    event_aligned_mua = AP_deconv_wfJF(event_aligned_mua, true);
end

% Wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, t_peri_event);

% Facecam movement
if exist('frame_movement', 'var')
    event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
        frame_movement(~isnan(facecam_t)), t_peri_event);
end

if task_dataset

    % Outcome (reward page 1, punish page 2)
    % (note incorrect outcome imprecise from signals, but looks good)
    event_aligned_outcome = zeros(size(t_peri_event, 1), size(t_peri_event, 2), 2);

    event_aligned_outcome(trial_outcome == 1, :, 1) = ...
        (cell2mat(arrayfun(@(x) ...
        histcounts(reward_t_timeline, t_peri_event_bins(x, :)), ...
        find(trial_outcome == 1), 'uni', false))) > 0;

    event_aligned_outcome(trial_outcome == -1, :, 2) = ...
        (cell2mat(arrayfun(@(x) ...
        histcounts(signals_events.responseTimes, t_peri_event_bins(x, :)), ...
        find(trial_outcome == -1), 'uni', false))) > 0;

end

%% Regress cortex to striatum and trial-align

if verbose;
    disp('Regressing cortex to striatum...');
end;

% Parameters for regression
% regression_params.use_svs = 1:100;
% regression_params.skip_seconds = 20;
% regression_params.upsample_factor = 1;
% regression_params.kernel_t = [-0.5, 0.5];
% regression_params.zs = [false, false];
% regression_params.cvfold = 5;
% regression_params.use_constant = true;
% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
% Get time points to bin
sample_rate = framerate * regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds, 1)):1 / sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds, 1, 'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins) / 2;

% Resample deconvolved fluorescence
fVdf_deconv_resample = interp1(frame_t, fVdf_deconv', time_bin_centers)';

% Bin spikes to match widefield frames
binned_spikes = nan(numel(theseNeurons), length(time_bin_centers));
binned_spikes_std = nan(numel(theseNeurons), length(time_bin_centers));

for curr_depth = 1:numel(theseNeurons)
    thisNeuron = theseNeurons(curr_depth);
    curr_spike_times = spike_times_timeline(spike_templates == thisNeuron);

    % (skip if no spikes at this depth)
    if isempty(curr_spike_times)
        continue
    end

    binned_spikes(curr_depth, :) = histcounts(curr_spike_times, time_bins);

    if ~isnan(nanstd(binned_spikes(curr_depth, :), [], 2)) && nanstd(binned_spikes(curr_depth, :), [], 2) ~= 0
        binned_spikes_std(curr_depth, :) = binned_spikes(curr_depth, :) ./ ...
            nanstd(binned_spikes(curr_depth, :), [], 2);
    else
        binned_spikes_std(curr_depth, :) = binned_spikes(curr_depth, :);
    end
end
%binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

% (filter MUA if selected)
if filter_mua
    binned_spikes = AP_deconv_wfJF(binned_spikes, true);
    binned_spikes_std = binned_spikes ./ nanstd(binned_spikes, [], 2);
end

% Load lambda from previously estimated and saved
lambda_fn = 'E:\for_Julie\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal, {ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day, ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

% Regress cortex to striatum
lambda = 200;
[ctx_str_k, ctxpred_spikes_std, explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs, :), ...
    binned_spikes_std, kernel_frames, lambda, ...
    regression_params.zs, regression_params.cvfold, ...
    true, regression_params.use_constant);

% Recast the k's into the master U
ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:, :, regression_params.use_svs), ...
    reshape(ctx_str_k{1}, size(ctx_str_k{1}, 1), []), U_master(:, :, regression_params.use_svs)), ...
    size(ctx_str_k{1}));

% Re-scale the prediction (subtract offset, multiply, add scaled offset)
ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})) .* ...
    nanstd(binned_spikes, [], 2) + ...
    nanstd(binned_spikes, [], 2) .* squeeze(ctx_str_k{end});

event_aligned_mua_ctxpred = ...
    interp1(time_bin_centers, ctxpred_spikes', t_peri_event) ./ raster_sample_rate;

%% Regress cortex to wheel velocity and speed and trial-align

if verbose;
    disp('Regressing cortex to wheel...');
end;

wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps, wheel_velocity, time_bin_centers);
wheel_velspeed_resample = [wheel_velocity_resample; abs(wheel_velocity_resample)];
wheel_velspeed_resample_std = wheel_velspeed_resample ./ std(wheel_velocity_resample);

[ctx_wheel_k, predicted_wheel_velspeed_std, explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs, :), ...
    wheel_velspeed_resample_std, kernel_frames, lambda, ...
    regression_params.zs, regression_params.cvfold, ...
    false, false);

predicted_wheel_velspeed = predicted_wheel_velspeed_std .* ...
    std(wheel_velocity_resample);

% Recast the k's into the master U
ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:, :, regression_params.use_svs), ...
    reshape(ctx_wheel_k, size(ctx_wheel_k, 1), []), U_master(:, :, regression_params.use_svs)), ...
    size(ctx_wheel_k));

event_aligned_wheel_ctxpred = ...
    interp1(time_bin_centers, predicted_wheel_velspeed(1, :)', t_peri_event);

%% Regress task to cortex/striatum/cortex-predicted striatum

if task_dataset

    if verbose
        disp('Regressing task to neural data...');
    end

    % Build regressors (only a subset of these are used)

    % Stim regressors
    unique_stim = unique(contrasts(contrasts > 0).*sides');
    stim_contrastsides = ...
        signals_events.trialSideValues(1:length(stimOn_times))' .* ...
        signals_events.trialContrastValues(1:length(stimOn_times))';

    stim_regressors = zeros(length(unique_stim), length(time_bin_centers));
    for curr_stim = 1:length(unique_stim)
        curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
        stim_regressors(curr_stim, :) = histcounts(curr_stim_times, time_bins);
    end

    % Stim move regressors (one for each stim when it starts to move)
    stim_move_regressors = zeros(length(unique_stim), length(time_bin_centers));
    for curr_stim = 1:length(unique_stim)

        % (find the first photodiode flip after the stim azimuth has
        % moved past a threshold)

        curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
            stim_contrastsides == unique_stim(curr_stim));

        azimuth_move_threshold = 5; % degrees to consider stim moved
        stim_move_times_signals = ...
            signals_events.stimAzimuthTimes( ...
            abs(signals_events.stimAzimuthValues-90) > azimuth_move_threshold);
        curr_stim_move_times_signals = arrayfun(@(x) ...
            stim_move_times_signals(find(stim_move_times_signals > ...
            curr_stimOn_times(x), 1)), 1:length(curr_stimOn_times));

        curr_stim_move_times_photodiode = arrayfun(@(x) ...
            photodiode_flip_times(find(photodiode_flip_times > ...
            curr_stim_move_times_signals(x), 1)), 1:length(curr_stim_move_times_signals));

        stim_move_regressors(curr_stim, :) = histcounts(curr_stim_move_times_photodiode, time_bins);

    end

    % Stim center regressors (one for each stim when it's stopped during reward)
    unique_contrasts = unique(contrasts(contrasts > 0));

    stim_center_regressors = zeros(length(unique_contrasts), length(time_bin_centers));
    for curr_contrast = 1:length(unique_contrasts)

        % (find the last photodiode flip before the reward)
        curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
            abs(stim_contrastsides) == unique_contrasts(curr_contrast));

        curr_reward_times = arrayfun(@(x) ...
            reward_t_timeline(find(reward_t_timeline > ...
            curr_stimOn_times(x), 1)), 1:length(curr_stimOn_times));

        curr_prereward_photodiode_times = arrayfun(@(x) ...
            photodiode_flip_times(find(photodiode_flip_times < ...
            curr_reward_times(x), 1, 'last')), 1:length(curr_reward_times));

        stim_center_regressors(curr_contrast, :) = histcounts(curr_prereward_photodiode_times, time_bins);

    end
%if task_dataset
    % Move onset regressors (L/R)
    move_onset_regressors = zeros(2, length(time_bin_centers));
    move_onset_regressors(1, :) = histcounts(wheel_move_time(trial_choice == -1), time_bins);
    move_onset_regressors(2, :) = histcounts(wheel_move_time(trial_choice == 1), time_bins);

    % Go cue regressors - separate for early/late move
    % (using signals timing - not precise but looks good)
    % (for go cue only on late move trials)
    %         go_cue_regressors = histcounts( ...
    %             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
    % (for go cue with early/late move trials)
    go_cue_regressors = zeros(1, length(time_bin_centers));
    go_cue_regressors(1, :) = histcounts( ...
        signals_events.interactiveOnTimes(stim_to_move <= 0.5), time_bins);
    go_cue_regressors(2, :) = histcounts( ...
        signals_events.interactiveOnTimes(stim_to_move > 0.5), time_bins);

    % Outcome regressors
    % (using signals timing - not precise but looks good)
    % (regressors for hit only)
    %         outcome_regressors = histcounts(reward_t_timeline,time_bins);
    % (regressors for both hit and miss)
    outcome_regressors = zeros(2, length(time_bin_centers));
    outcome_regressors(1, :) = histcounts( ...
        reward_t_timeline, time_bins);
    outcome_regressors(2, :) = histcounts( ...
        signals_events.responseTimes(trial_outcome == -1), time_bins);

    % Concatenate selected regressors, set parameters
    task_regressors = {stim_regressors; move_onset_regressors; go_cue_regressors; outcome_regressors};
    task_regressor_labels = {'Stim', 'Move onset', 'Go cue', 'Outcome'};

    task_t_shifts = { ...
        [0, 0.5]; ... % stim
        [-0.5, 1]; ... % move
        [0, 0.5]; ... % go cue
        [0, 0.5]}; % outcome
% else
%      task_regressors = {stim_regressors};
%     task_regressor_labels = {'Stim'};
% 
%     task_t_shifts = { ...
%         [0, 0.5]}; % stim
% end
    task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
        round(x(2)*(sample_rate)), task_t_shifts, 'uni', false);
    lambda = 0;
    zs = [false, false];
    cvfold = 5;
    use_constant = false;
    return_constant = false;

    % Regression task -> MUA
    baseline = nanmean(reshape(event_aligned_mua(:, t < 0, :), [], ...
        size(event_aligned_mua, 3))*raster_sample_rate, 1)';
    activity = single(binned_spikes) - baseline;

    [mua_taskpred_k, mua_taskpred_long, mua_taskpred_expl_var, mua_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors, activity, task_regressor_sample_shifts, ...
        lambda, zs, cvfold, return_constant, use_constant);

    mua_taskpred = ...
        interp1(time_bin_centers, mua_taskpred_long', t_peri_event) ./ raster_sample_rate;

    mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers, mua_taskpred_reduced_long(:, :, x)', ...
        t_peri_event)./raster_sample_rate, permute(1:length(task_regressors), [1, 3, 4, 2]), 'uni', false));

    % Regression task -> MUA-ctxpred
    baseline = nanmean(reshape(event_aligned_mua_ctxpred(:, t < 0, :), [], ...
        size(event_aligned_mua_ctxpred, 3))*raster_sample_rate, 1)';
    activity = single(ctxpred_spikes) - baseline;

    [mua_ctxpred_taskpred_k, mua_ctxpred_taskpred_long, mua_ctxpred_taskpred_expl_var, mua_ctxpred_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors, activity, task_regressor_sample_shifts, ...
        lambda, zs, cvfold, return_constant, use_constant);

    mua_ctxpred_taskpred = ...
        interp1(time_bin_centers, mua_ctxpred_taskpred_long', t_peri_event) ./ raster_sample_rate;

    mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers, mua_ctxpred_taskpred_reduced_long(:, :, x)', ...
        t_peri_event)./raster_sample_rate, permute(1:length(task_regressors), [1, 3, 4, 2]), 'uni', false));

    % Regression task -> (master U, deconvolved) fluor
    event_aligned_V_deconv = AP_deconv_wfJF(event_aligned_V);
    fVdf_deconv_resample_recast = ChangeU(Udf_aligned, fVdf_deconv_resample, U_master);

    baseline = nanmean(reshape(event_aligned_V_deconv(:, t < 0, :), [], size(event_aligned_V_deconv, 3)))';
    activity = single(fVdf_deconv_resample_recast(use_components, :)) - baseline;

    [fluor_taskpred_k, fluor_taskpred_long, fluor_taskpred_expl_var, fluor_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors, activity, task_regressor_sample_shifts, ...
        lambda, zs, cvfold, return_constant, use_constant);

    fluor_taskpred = ...
        interp1(time_bin_centers, fluor_taskpred_long', t_peri_event);

    fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers, fluor_taskpred_reduced_long(:, :, x)', ...
        t_peri_event), permute(1:length(task_regressors), [1, 3, 4, 2]), 'uni', false));

end

%% Store everything into structure



trial_data.trial_info_all = trial_info;

trial_data.fluor_all = event_aligned_V(use_trials, :, :, :);
trial_data.mua_all = event_aligned_mua(use_trials, :, :, :);
% trial_data.mua_move_all = event_aligned_move_mua(use_trials, :, :, :);
% trial_data.mua_outcome_all = event_aligned_outcome_mua(use_trials, :, :, :);

trial_data.ctx_str_k_all = ctx_str_k_recast;
trial_data.mua_ctxpred_all = event_aligned_mua_ctxpred(use_trials, :, :, :);

trial_data.wheel_all = event_aligned_wheel(use_trials, :, :);
if exist('event_aligned_movement', 'var')
    trial_data.movement_all = event_aligned_movement(use_trials, :, :);
end

trial_data.ctx_wheel_k_all = ctx_wheel_k_recast;
trial_data.wheel_ctxpred_all = event_aligned_wheel_ctxpred(use_trials, :, :);
   trial_data.allGroups = allGroups;
    trial_data.allDepths = allDepths;
    trial_data.allAnimals = allAnimals;
    trial_data.allDays = allDays;
    trial_data.allExperiments = allExperiments;
    trial_data.allNeurons = theseNeurons;
    
if task_dataset
 
    
    trial_data.outcome_all = event_aligned_outcome(use_trials, :, :);

    trial_data.mua_taskpred_k_all = mua_taskpred_k;
    trial_data.mua_taskpred_all = mua_taskpred(use_trials, :, :, :);
    trial_data.mua_taskpred_reduced_all = mua_taskpred_reduced(use_trials, :, :, :);
    trial_data.mua_taskpred_expl_var_total_all = mua_taskpred_expl_var.total;
    trial_data.mua_taskpred_expl_var_partial_all = mua_taskpred_expl_var.partial;

    trial_data.mua_ctxpred_taskpred_k_all = mua_ctxpred_taskpred_k;
    trial_data.mua_ctxpred_taskpred_all = mua_ctxpred_taskpred(use_trials, :, :, :);
    trial_data.mua_ctxpred_taskpred_reduced_all = mua_ctxpred_taskpred_reduced(use_trials, :, :, :);
    trial_data.mua_ctxpred_taskpred_expl_var_total_all = mua_ctxpred_taskpred_expl_var.total;
    trial_data.mua_ctxpred_taskpred_expl_var_partial_all = mua_ctxpred_taskpred_expl_var.partial;

    trial_data.fluor_taskpred_k_all = fluor_taskpred_k;
    trial_data.fluor_taskpred_all = fluor_taskpred(use_trials, :, :, :);
    trial_data.fluor_taskpred_reduced_all = fluor_taskpred_reduced(use_trials, :, :, :);
    trial_data.fluor_taskpred_expl_var_total_all = fluor_taskpred_expl_var.total;
    trial_data.fluor_taskpred_expl_var_partial_all = fluor_taskpred_expl_var.partial;

end

if verbose
    disp('Done getting trial data.');
end
