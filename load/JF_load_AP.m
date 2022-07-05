
%% Find experiments with the task + cortical widefield + striatal ephys
close all;
myPaths;

animals={'JF070'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'oice'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

% allen_atlas_path = [allenAtlasPath, filesep, 'allenCCF/'];
% st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
% curr_plot_structure = find(strcmp(st.acronym, 'GPe'));
% histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
% load(histoFile)
% probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
% load(probe2ephysFile)
% min(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))
% max(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))

%% Load data from experiment 

curr_day = 6; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
date = thisDay;
verbose = false; % display load progress and some info figures
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true;

site = 	1;%1,1; 2,4; 3,7
recording = []; 
experiment = 2;
loadClusters = 0;
[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_ap',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = contains(ephysAPfile, '_g');
if isSpikeGlx
     [ephysKSfile,~] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys',site,recording);
    if isempty(dir([ephysKSfile filesep 'sync.mat']))
        sync = syncFT(ephysAPfile, 385, ephysKSfile);
    end
end

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
if ~isempty(qMetricsExist)

load(fullfile(savePath, 'qMetric.mat'))
load(fullfile(savePath, 'param.mat'))
bc_getQualityUnitType;
end
clearvars unitType 
AP_load_experimentJF;
curr_shank=NaN;

trial_conditions(~ismember(trial_conditions(:,1), [4,6,7]),1) = 0;
AP_cellrasterJF({stimOn_times(~isnan(stimOn_times)), stimOn_times(~isnan(stimOn_times)), stimOn_times(~isnan(stimOn_times))}, ...
    {trial_conditions(~isnan(stimOn_times),1), trial_conditions(~isnan(stimOn_times),2),(trial_conditions(~isnan(stimOn_times),2)+1).*(trial_conditions(~isnan(stimOn_times),1)+1)});

AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))'}, ...
    {trial_conditions(:,1),trial_conditions(:,2), ...
    trial_conditions(:,3)});

AP_cellrasterJF({stimOn_times(~isnan(stimOn_times)), stimOn_times(~isnan(stimOn_times))}, ...
    {trial_conditions(~isnan(stimOn_times),1), trial_conditions(~isnan(stimOn_times),2)});

trial_conditions(~ismember(trial_conditions(:,1), [4,6,7]),1) = 0;
AP_cellrasterJF({stimOn_times(~isnan(stimOn_times)), stimOn_times(~isnan(stimOn_times))}, ...
    {trial_conditions(~isnan(stimOn_times),1), trial_conditions(~isnan(stimOn_times),2)});

AP_cellrasterJF({stimOn_times(ismember(trial_conditions(:,1), [4,6,7])), stimOn_times(ismember(trial_conditions(:,1), [4,6,7]))}, ...
    {trial_conditions(ismember(trial_conditions(:,1), [4,6,7]),1), trial_conditions(ismember(trial_conditions(:,1), [4,6,7]),2)});


AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))',stimOn_times}, ...
    {trial_conditions(:,1),trial_conditions(:,2), ...
    trial_conditions(:,3),movement_after200ms_and_type});

AP_cellrasterJF({stimOn_times, stimOn_times(ismember(trial_conditions(:,1), [4,6,7]))}, ...
    {trial_conditions(:,2), trial_conditions(ismember(trial_conditions(:,1), [4,6,7]),1)});




%% Open widefield viewer

% This is to get a feel for the widefield data, it shows the mouse cameras
% (upper plots), the widefield fluorescence (lower left), and stimulus on
% lower right (darker color = higher contrast, border color = status for
% immobile/interactive/rewarded/punished)

fVdf_deconv = AP_deconv_wf(fVdf);
AP_wfmovies_choiceworld(Udf,fVdf_deconv,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t,signals_events)


%% Regress cortex to striatum units

% (set the sample rate - match to framerate)
sample_rate = (1/median(diff(frame_t)));

% (get time bins to match fluorescence and spikes)
% (skip the first/last n seconds to do this, sometimes imaging artifacts)
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (only use spikes/units within time bins and in striatum)
use_spikes = spike_times_timeline >= time_bins(1) & spike_times_timeline <= time_bins(end) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
template_spike_n = accumarray(spike_templates(use_spikes),1,[size(templates,1),1]);
use_templates = find(template_spike_n > 0);

% Bin spikes for each unit
binned_spikes = zeros(length(use_templates),length(time_bins)-1);
for curr_template_idx = 1:length(use_templates)
    curr_template = use_templates(curr_template_idx);
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);
    binned_spikes(curr_template_idx,:) = histcounts(curr_spike_times,time_bins);
end

% Deconvolve fluorescence and match to spike times
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

% Set parameters for regression
use_svs = 1:50;
kernel_t = [-0.2,0.2];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 10;
zs = [false,false];
cvfold = 5;

% Regress cortex to units (using deconvolved fluorescence)
% (std units to have common lambda regularization)
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(use_svs,:), ...
    binned_spikes_std, ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
k_px = arrayfun(@(x) svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,x)),1:size(k,3),'uni',false);
% Concatenate and flip time (to have fluroescence lead->lag spikes)
k_px = flip(cat(4,k_px{:}),3);

% Plot cortex -> striatum unit kernels
% (this uses a 4D matrix viewer: scroll wheel changes dim 3 (time here),
% up/down arrows changes dim 4 (units here))
t_label = cellfun(@(x) ['Fluor relative to spike: ' num2str(x) ' s'], ...
    num2cell(kernel_frames/sample_rate),'uni',false);
AP_image_scroll(k_px,t_label)
axis image;
caxis([-0.005,0.005])
colormap(brewermap([],'*RdBu'));

%% Regress task events to ephys units
% messy at the moment - I can explain more if you want

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get event-aligned activity
raster_window = [-0.5,2];
upsample_factor = 1;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;

t_peri_event = bsxfun(@plus,use_align,t);
t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

%%% Trial-align wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

%%% Trial-align facecam movement
event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
    frame_movement(~isnan(facecam_t)),t_peri_event);

%%% Trial-align outcome (reward page 1, punish page 2)
% (note incorrect outcome imprecise from signals, but looks good)
event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);

event_aligned_outcome(trial_outcome == 1,:,1) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(reward_t_timeline,t_bins(x,:)), ...
    find(trial_outcome == 1)','uni',false))) > 0;

event_aligned_outcome(trial_outcome == -1,:,2) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(signals_events.responseTimes,t_bins(x,:)), ...
    find(trial_outcome == -1)','uni',false))) > 0;

% Pick trials to keep
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = signals_events.trialSideValues(1:n_trials) == -1;
R_trials = signals_events.trialSideValues(1:n_trials) == 1;

D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);

D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
D.repeatNum = ones(sum(use_trials),1);

D.outcome = reshape(trial_outcome(use_trials),[],1);

%%% Regress task to cortex/striatum/cortex-predicted striatum

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get reaction time for building regressors
[move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
move_idx(~move_trial) = NaN;
move_t = nan(size(move_idx));
move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';

% Stim regressors
unique_stim = unique(contrasts(contrasts > 0).*sides');
stim_contrastsides = ...
    signals_events.trialSideValues(1:length(stimOn_times)).* ...
    signals_events.trialContrastValues(1:length(stimOn_times));

stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
    stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
end

% Stim move regressors (one for each stim when it starts to move)
stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    
    % (find the first photodiode flip after the stim azimuth has
    % moved past a threshold)
    
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
        stim_contrastsides == unique_stim(curr_stim));
    
    azimuth_move_threshold = 5; % degrees to consider stim moved
    stim_move_times_signals = ...
        signals_events.stimAzimuthTimes( ...
        abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
    curr_stim_move_times_signals = arrayfun(@(x) ...
        stim_move_times_signals(find(stim_move_times_signals > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_stim_move_times_photodiode = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times > ...
        curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
    
    stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
    
end

% Stim center regressors (one for each stim when it's stopped during reward)
unique_contrasts = unique(contrasts(contrasts > 0));
stim_contrasts = ...
    signals_events.trialContrastValues(1:length(stimOn_times));

stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
for curr_contrast = 1:length(unique_contrasts)
    
    % (find the last photodiode flip before the reward)
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
        stim_contrasts == unique_contrasts(curr_contrast));
    
    curr_reward_times = arrayfun(@(x) ...
        reward_t_timeline(find(reward_t_timeline > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_prereward_photodiode_times = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times < ...
        curr_reward_times(x),1,'last')),1:length(curr_reward_times));
    
    stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
    
end

% Move onset regressors (L/R)
move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));

move_onset_regressors = zeros(2,length(time_bin_centers));
move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);

% Move onset x stim regressors (one for each contrast/side)
move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
    arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
    1:length(unique_stim),'uni',false);

move_onset_stim_regressors = zeros(10,length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    move_onset_stim_regressors(curr_stim,:) = ...
        histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
end

% Move ongoing regressors (L/R choice for duration of movement)
wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);

move_stopped_t = 0.5;
move_stopped_samples = round(sample_rate*move_stopped_t);
wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
    ones(1,move_stopped_samples),'full') > 0;
wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);

move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_L_absolute','uni',false));
move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_R_absolute','uni',false));

move_ongoing_regressors = zeros(2,length(time_bin_centers));
move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
move_ongoing_regressors(2,move_ongoing_R_samples) = 1;

% Go cue regressors - separate for early/late move
% (using signals timing - not precise but looks good)
if length(signals_events.interactiveOnTimes) ~= length(move_t)
    error('Different number of interactive ons and move times')
end

go_cue_regressors = zeros(2,length(time_bin_centers));
go_cue_regressors(1,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
go_cue_regressors(2,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t > 0.5),time_bins);

% Outcome regressors
% (using signals timing - not precise but looks good)
outcome_regressors = zeros(2,length(time_bin_centers));

outcome_regressors(1,:) = histcounts( ...
    reward_t_timeline,time_bins);
outcome_regressors(2,:) = histcounts( ...
    signals_events.responseTimes(trial_outcome == -1),time_bins);

% Concatenate selected regressors, set parameters
regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;


%%%%%%% New: regression task -> striatum units

%%% Trial-align striatum units
event_aligned_unit = nan(length(stimOn_times),length(t),max(spike_templates));
t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
for curr_unit = 1:max(spike_templates)
    curr_spikes = spike_times_timeline(spike_templates == curr_unit);
    event_aligned_unit(:,:,curr_unit) = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes,t_bins(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
    AP_print_progress_fraction(curr_unit,max(spike_templates));
end

%%% Binned unit activity across the experiment
binned_spikes = zeros(max(spike_templates),length(time_bin_centers));
for curr_unit = 1:max(spike_templates)
    curr_spike_times = spike_times_timeline(spike_templates == curr_unit);
    binned_spikes(curr_unit,:) = histcounts(curr_spike_times,time_bins);
    AP_print_progress_fraction(curr_unit,max(spike_templates));
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

%%% Task > striatal unit regression
baseline = nanmean(reshape(event_aligned_unit(:,t < 0,:),[], ...
    size(event_aligned_unit,3))*raster_sample_rate)';
binned_spikes_baselinesubtracted = binned_spikes - baseline;

[unit_taskpred_k,~,unit_expl_var,~] = ...
    AP_regresskernel(regressors,binned_spikes_baselinesubtracted, ...
    sample_shifts,lambda,zs,cvfold,return_constant,use_constant);

% Plot unit explained variance by depth
used_spikes = spike_times_timeline > time_bins(1) & ...
    spike_times_timeline < time_bins(end);
norm_spike_n = mat2gray(log(accumarray(spike_templates(used_spikes),1,[size(templates,1),1])+1));

use_partial_var = 2;
[~,max_regressor_idx] = max(unit_expl_var.partial(:,:,use_partial_var),[],2);
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

figure;

subplot(1,length(regressors)+2,1,'YDir','reverse'); hold on;
scatter3(norm_spike_n,template_depths, ...
    1:max(spike_templates),20, ...
    regressor_cols(max_regressor_idx,:),'filled')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title({animal,day,'Best regressor'})

subplot(1,length(regressors)+2,2,'YDir','reverse'); hold on;
curr_expl_var = unit_expl_var.total;
curr_expl_var_dotsize = 100*mat2gray(curr_expl_var,[0,1]) + 1;
curr_expl_var_dotsize(isnan(curr_expl_var) | curr_expl_var < -1 | curr_expl_var >= 1) = NaN;
scatter3(norm_spike_n,template_depths, ...
    1:max(spike_templates),curr_expl_var_dotsize, ...
    regressor_cols(max_regressor_idx,:),'filled')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title('Total expl var');

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressors)
    subplot(1,length(regressors)+2,curr_regressor+2,'YDir','reverse');
    hold on;
    curr_expl_var = unit_expl_var.partial(:,curr_regressor,use_partial_var);
    curr_nan = isnan(curr_expl_var) | curr_expl_var < 0 | curr_expl_var >= 1;
    curr_expl_var(curr_nan) = NaN;
    curr_expl_var = mat2gray(curr_expl_var);
    curr_expl_var(curr_nan) = NaN;
    scatter3(norm_spike_n,template_depths, ...
        1:max(spike_templates),curr_expl_var*50+1, ...
        regressor_cols(curr_regressor,:),'filled')
    title(regressor_labels{curr_regressor});
    xlabel('Normalized n spikes');
    ylabel('Depth (\mum)');
end













