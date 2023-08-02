
%% Find experiments with the task + cortical widefield + striatal ephys
%close all;
%clear all;
cl_myPaths;

animals={'JF093'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

% allen_atlas_path = [allenAtlasPath, filesep, 'allenCCF/'];
% st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
% curr_plot_structure = find(strcmp(st.acronym, 'GPe'));m
% histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
% load(histoFile)
% probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
% load(probe2ephysFile)
% min(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))
% max(probe_ccf(9).probe_depths(probe_ccf(9).trajectory_areas ==curr_plot_structure))

%% Load data from experiment 

curr_day = 1; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
date = thisDay;
verbose = false; % display load progress and some info figures
load_parts.cam=false;
load_parts.imaging=false;
load_parts.ephys=true;

site = 1;%1,1; 2,4; 3,7
recording = []; 
% keep experiment with max n trials (= most likely not aborted error or end
% shank mapping) QQ change this in the future 
n_trials = zeros(size(experiments(curr_day).experiment,2), 1);
for iExperiment = experiments(curr_day).experiment
    [block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, iExperiment, 'block');
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

experiment = 2;%experiments(curr_day).experiment(1);%find(n_trials == max(n_trials));
loadClusters = 0;
[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_includingCompressed',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
% if ~isempty(qMetricsExist)
% 
% [param, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(savePath);
% unitType = bc_getQualityUnitType(param, qMetric);
% end
%clearvars unitType 
%load_parts.cam = true;
%load_parts.cam = false;

JF_load_experiment;
curr_shank=NaN;
%AP_cellraster({stimOn_times, stimOn_times}, {trial_conditions(:,2), trial_conditions(:,3)})

trial_conditions_clean = trial_conditions;
trial_conditions_clean(trial_conditions(:,1)==4,1) = 100;%go1
trial_conditions_clean(trial_conditions(:,1)==12,1) = 101;%go2
trial_conditions_clean(trial_conditions(:,1)==6,1) = 102;%no go
trial_conditions_clean(trial_conditions(:,1)==11,1) = 103;%no like
trial_conditions_clean(trial_conditions(:,1)==13,1) = 104;%go like
theseImages = [100:102];
%trial_conditions(trial_conditions(:,1)>13,1) = trial_conditions(trial_conditions(:,1)>13,1)-13;

theseImages_trials = ismember(trial_conditions_clean(:,1), theseImages) & ismember(trial_conditions_clean(:,2), [-90]);
AP_cellrasterJF({stimOn_times(theseImages_trials), stimOn_times(theseImages_trials), stimOn_times(theseImages_trials)},...
    {trial_conditions_clean(theseImages_trials,1),...
    trial_conditions_clean(theseImages_trials,2), -89+trial_conditions_clean(theseImages_trials,1) + abs(trial_conditions_clean(theseImages_trials,2))})
%AP_cellraster({stimOn_times, stimOn_times}, {trial_conditions(:,2), trial_conditions(:,3)})
%AP_cellraster({stimOn_times}, {trial_conditions(:,1)})
%AP_cellraster({stimOn_times, stimOn_times}, {trial_conditions(:,1), trial_conditions(:,2)})
% 
% % passive, task stims 
% theseImages = [4,6,12];
% trial_conditions(trial_conditions(:,1)>13,1) = trial_conditions(trial_conditions(:,1)>13,1)-13;
% 
% theseImages_trials = ismember(trial_conditions(:,1), theseImages);
% 
% AP_cellrasterJF({stimOn_times(theseImages_trials & trial_conditions(:,2)==-90), stimOn_times(theseImages_trials & trial_conditions(:,2)==0),...
%     stimOn_times(theseImages_trials), stimOn_times(theseImages_trials)}, ...
%     {trial_conditions(theseImages_trials & trial_conditions(:,2)==-90,1), ...
%     trial_conditions(theseImages_trials & trial_conditions(:,2)==0,1),trial_conditions(theseImages_trials,2),...
% (trial_conditions(theseImages_trials,2)/-90)+(trial_conditions(theseImages_trials,1))});
% 
% 
% % go no go 
% correct_trials = trial_conditions(:,3)==1;
% correct_trials_no_repeat = trial_conditions(:,3)==1 & repeatOnIncorrect' == 1;
% AP_cellraster({stimOn_times, stimOn_times(correct_trials), stimOn_times(correct_trials),...
%     wheel_move_time,signals_events.responseTimes(n_trials)'}, ...
%     {trial_conditions(:,1), trial_conditions(correct_trials,1),...
%     trial_conditions(correct_trials,1), trial_conditions(:,2), ...
%     trial_conditions(:,3)});
% 
% % go go go 
% AP_cellraster({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials)', ...
%     [signals_events.responseTimes(crrectExpect)'; unexpected_reward_times']}, ...
%     {trial_conditions(:,1),trial_conditions(:,2), ...
%     trial_conditions(:,3), [trial_outcome_reward_omission'; 2*ones(size(unexpected_reward_times,2),1)]});
% % population raster 
% 
% experiment = 2;
% JF_load_experiment;
% curr_shank=NaN;
% theseImages = [4,6,12];
% trial_conditions(trial_conditions(:,1)>13,1) = trial_conditions(trial_conditions(:,1)>13,1)-13;
% theseImages_trials = ismember(trial_conditions(:,1), theseImages);
% AP_cellrasterJF({stimOn_times(theseImages_trials & trial_conditions(:,2)==-90), stimOn_times(theseImages_trials & trial_conditions(:,2)==0),...
%     stimOn_times(theseImages_trials), stimOn_times(theseImages_trials)}, ...
%     {trial_conditions(theseImages_trials& trial_conditions(:,2)==-90,1), ...
%     trial_conditions(theseImages_trials& trial_conditions(:,2)==0,1),trial_conditions(theseImages_trials,2),...
% (trial_conditions(theseImages_trials,2)/-90)+(trial_conditions(theseImages_trials,1))});
% % 
% % % unexpected_reward_times 
% % %   trial_conditions = ...
% % %                 [signals_events.stimulusTypeValues(n_trials(1):n_trials(end))', ...
% % %                 trial_choice(n_trials(1):n_trials(end))', trial_outcome(n_trials(1):n_trials(end))',...
% % %                 trial_outcome_expected(n_trials(1):n_trials(end)), trial_outcome_reward_omission(n_trials(1):n_trials(end))'];
% %        
% % 
% % % get behav measures, get pupil 
% % [instHit_rate, instCR_rate ] = JF_getBehavArousalMeasures(stimIDs, trial_choice, stimOn_times);
% % % instHit_rate = movmean(double(~isnan(stim_to_move(ismember(stimIDs(1:length(n_trials)), [1, 2])))), 20);
% % % instCR_rate = movmean(double(isnan(stim_to_move(ismember(stimIDs(1:length(n_trials)), [3])))), 20);
% % 
% % GoTrials = trial_conditions(:,1) ==1;
% % JF_cellraster({stimOn_times(GoTrials),wheel_move_time,signals_events.responseTimes(n_trials)'}, ...
% %     {trial_conditions(GoTrials,2),trial_conditions(:,2), ...
% %     trial_conditions(:,3)}, {stim_to_move, stim_to_feedback-stim_to_move, -stim_to_feedback+stim_to_move},...
% %     [],{instHit_rate, instCR_rate});
% % 
% % noGoTrials = trial_conditions(:,1) ==3;
% % JF_cellraster({stimOn_times(noGoTrials),wheel_move_time,signals_events.responseTimes(n_trials)'}, ...
% %     {trial_conditions(noGoTrials,1),trial_conditions(:,2), ...
% %     trial_conditions(:,3)}, {stim_to_move, stim_to_feedback-stim_to_move, -stim_to_feedback+stim_to_move},...
% %     [],{instHit_rate, instCR_rate});
% % % 
% % % % % go no go passive
% % % % 
% % % %  trial_conditions(ismember(trial_conditions(:,1), [4]),1) = 4; % go 1
% % % % 
% % % % 
% % % %  trial_conditions(ismember(trial_conditions(:,1), [7]),1) = 7; % go 2
% % % %  trial_conditions(ismember(trial_conditions(:,1), [10]),1) = 1; % no go
% % % %  trial_conditions(ismember(trial_conditions(:,1), [6]),1) = 10; % no go
% % % %   trial_conditions(~ismember(trial_conditions(:,1), [4,7,10]),1) = 1;
% % % % %thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions, [4,7,10]) & trial_conditions(:,2)~=90;
% % % % thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions(:,2), [-90,0]) & ismember(trial_conditions(:,1), [4,7,10]);
% % % % 
% % % % 
% % AP_cellrasterJF({stimOn_times, stimOn_times, stimOn_times}, ...
% %     {trial_conditions(:,1), trial_conditions(:,2),...
% % (trial_conditions(:,2)/-90)+(trial_conditions(:,1))});
% % % 
% % %  trial_conditions(ismember(trial_conditions(:,1), [3]),1) = 4; % go 1
% % %  trial_conditions(ismember(trial_conditions(:,1), [4]),1) = 4; % go 1
% % % 
% % % 
% % %  trial_conditions(ismember(trial_conditions(:,1), [7]),1) = 7; % go 2
% % %  trial_conditions(ismember(trial_conditions(:,1), [8]),1) = 7; % go 2
% % %  trial_conditions(ismember(trial_conditions(:,1), [10]),1) = 1; % no go
% % %  trial_conditions(ismember(trial_conditions(:,1), [6]),1) = 10; % no go
% % %  trial_conditions(ismember(trial_conditions(:,1), [5]),1) = 10; % no go
% % %  trial_conditions(~ismember(trial_conditions(:,1), [4,7,10]),1) = 1;
% % % %thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions, [4,7,10]) & trial_conditions(:,2)~=90;
% % % thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions(:,2), [-90,0]) & ismember(trial_conditions(:,1), [4,7,10]);
% % % 
% % % %  trial_conditions(ismember(trial_conditions(:,1), [1,2,3,4,18,19,20,21,22]),1) = 4; % go 1
% % % %  trial_conditions(ismember(trial_conditions(:,1), [5,6,10,11,12,13,14,15,16,17]),1) = 10; % go 1
% % % %   trial_conditions(ismember(trial_conditions(:,1), [7,8,9]),1) = 7; % go 1
% % % % thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions(:,2), [-90,0]) & ismember(trial_conditions(:,1), [4,7,10]);
% % % % 
% % % 
% % % AP_cellrasterJF({stimOn_times(thisIndex), stimOn_times(thisIndex), stimOn_times(thisIndex)}, ...
% % %     {trial_conditions(thisIndex,1), trial_conditions(thisIndex,2),...
% % % (trial_conditions(thisIndex,2)/-90)+(trial_conditions(thisIndex,1))});
% % % 
% % % thisIndex = ~isnan(stimOn_times(1:size(trial_conditions,1))) & ismember(trial_conditions(:,1), [4,6,7]);
% % % 
% % % AP_cellrasterJF({stimOn_times(thisIndex), stimOn_times(thisIndex), stimOn_times(thisIndex)}, ...
% % %     {trial_conditions(thisIndex,1), trial_conditions(thisIndex,2),...
% % % (trial_conditions(thisIndex,2)/-90)+(trial_conditions(thisIndex,1))});
% % % 
% % % % imageworld passive
% % % trial_conditions(ismember(trial_conditions(:,1), [1:3:66]),2) = -90;
% % % trial_conditions(ismember(trial_conditions(:,1), [2:3:66]),2) = 0;
% % % trial_conditions(ismember(trial_conditions(:,1), [3:3:66]),2) = 90;
% % % trial_conditions(ismember(trial_conditions(:,1), [4,26,48])) = 4;
% % % trial_conditions(ismember(trial_conditions(:,1), [7,29,51])) = 7;
% % % trial_conditions(ismember(trial_conditions(:,1), [10,32,54])) = 10;
% % % thisIndex = ismember(trial_conditions(:,1), [4,7,10]) & trial_conditions(:,2)~=-90; %ismember(trial_conditions(:,1), [4,7,10]);
% % % 
% % % AP_cellrasterJF({stimOn_times(thisIndex), stimOn_times(thisIndex), stimOn_times(thisIndex)}, ...
% % %     {trial_conditions(thisIndex,1), trial_conditions(thisIndex,2),...
% % % (trial_conditions(thisIndex,2)/-90)+(trial_conditions(thisIndex,1))});
% % % 
% % % 
% % % AP_cellrasterJF({stimOn_times}, {trial_conditions(:,1)})
% % % % task 
% % AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(1:end-1)'}, ...
% %     {trial_conditions(:,1),trial_conditions(:,2), ...
% %     trial_conditions(:,3)});
% % 
% % % % 
% % % % AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))',stimOn_times}, ...
% % % %     {trial_conditions(:,1),trial_conditions(:,2), ...
% % % %     trial_conditions(:,3), movement_after200ms_and_type});
% % % 
% % % 
% % % % PSTH more in depth 
% % % figure(1)
% % % thisTemplate = 79;
% % % raster_window = [-0.5, 2];
% % % align_times = stimOn_times(stimIDs==3);
% % % align_group = [];
% % % color_by = trial_conditions(stimIDs==3,3);
% % % psth_bin_size = 0.001;
% % % sort_by = [];%stim_to_move(stimIDs==3);
% % % plot_me = true;
% % % [curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
% % %     thisTemplate, raster_window, psth_bin_size, align_times, align_group,...
% % %    sort_by, color_by, plot_me);
% % % title(['unit' num2str(thisTemplate)])
% % % 
% % % figure(2)
% % % 
% % % raster_window = [-0.5, 2];
% % % align_times = stimOn_times(stimIDs==2);
% % % align_group = [];
% % % color_by = trial_conditions(stimIDs==2,3);
% % % psth_bin_size = 0.001;
% % % sort_by = [];%stim_to_move(stimIDs==3);
% % % plot_me = true;
% % % [curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
% % %     thisTemplate, raster_window, psth_bin_size, align_times, align_group,...
% % %    sort_by, color_by, plot_me);
% % % title(['unit' num2str(thisTemplate)])


%% example cells 
% CP: 110, 1, 1 units 241, 185, 235, 223, 215, 325 
stimC = [224, 85, 159]./256;
unique_templates = unique(spike_templates);
%responds 
thisUnit = 46
colorMtx =  bc_colors(3, 'w');
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
align_times = stimOn_times;
[curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
    thisUnit , raster_window, psth_bin_size, align_times, []);
               figure();
subplot(3,1,[1,2])
s = scatter(t(raster_x(1:1:end)), raster_y(1:1:end), 3, [0,0,0], 'filled', 'MarkerFaceAlpha', 0.1, ...
    'MarkerEdgeAlpha', 0.1);
s.MarkerFaceAlpha = 0.2;
s.AlphaData = 0.1* ones(size(raster_y,1),1);
s.MarkerFaceAlpha = 'flat';
ylabel('trial #')
xticklabels({''})
makepretty;

ylim([1 1080])
yLim = ylim; 
line([0, 0], [yLim(1), yLim(2)], 'Color', stimC, 'LineWidth', 2)
xlim([-0.15, 0.3])

psth_bin_size = 0.01;
align_times = stimOn_times;
[curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
    thisUnit , raster_window, psth_bin_size, align_times, []);
subplot(3,1,[3])
plot(t, curr_psth .* 1/psth_bin_size,'Color', [0,0,0])
xlabel('time from stim onset (s)')
ylabel('spikes/s')
makepretty;
hold on;
yLim = ylim; 
line([0, 0], [yLim(1), yLim(2)], 'Color', stimC, 'LineWidth', 2)
xlim([-0.15, 0.3])
% selectivity
