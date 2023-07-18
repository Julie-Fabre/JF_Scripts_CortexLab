
%example session

% extract TANs, FSIs, MSNs


% plot activity in first x trails, and below add lapses / gos (two traces)
% add pupil too ? 

%DMS: 93 probe 1 [3,1], *2 [1,1]*, (3, bott)[2,1], 7(top)[7,1], *8[8,1]*,
%9(top)[8,2], 10(top)[9,1]
%GPe: 93 probe 4[6,1], *5[5,1]*, *6[4,1]*, 7[7,1], 9[8,2]
%GPi: 93 probe 10[9,1]
%SNr: 93 probe *14[12,2]*


% load example dataset 
cl_myPaths;
animal='JF093';
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
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

experiment = 1;%experiments(curr_day).experiment(1);%find(n_trials == max(n_trials));
loadClusters = 0;
[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_includingCompressed',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end

ephysDirPath = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site, recording);

JF_load_experiment;

% get good units 
savePathQMetric = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePathQMetric, '*qMetric*'));
if ~isempty(qMetricsExist)
    [param, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(savePathQMetric);
    unitType = bc_getQualityUnitType(param, qMetric);
end

% get unit types 
savePathEProp= fullfile(ephysDirPath, 'ephysProperties');
ephysPropExist = dir(fullfile(savePathEProp, '*ephysProperties*'));
if ~isempty(ephysPropExist)
    [paramEP, ephysProperties, acg] = bc_loadSavedProperties(savePathEProp);
    figure();
    scatter(ephysProperties.templateDuration, ephysProperties.postSpikeSuppression)
    msn = ephysProperties.templateDuration >= 400 & ephysProperties.postSpikeSuppression < 40;
    fsi = ephysProperties.templateDuration < 400;
    tan = ephysProperties.templateDuration >= 400 & ephysProperties.postSpikeSuppression >= 40;
end
% get units in striatum 
depthStriatum = 1500; 
goodStriatalUnits = (unitType == 1) & template_depths > 1500;

% raster msns, fsis, tans for 20 trials 
allGoodMSNs = find(goodStriatalUnits & msn);
allGoodTANs = find(goodStriatalUnits & tan);
allGoodFSIs = find(goodStriatalUnits & fsi);
allGoodCells = [allGoodMSNs; allGoodTANs; allGoodFSIs];
iChunk = 2;
chunkSize = 10;
trialStart = stimOn_times(1+(iChunk-1)*chunkSize);
trialStop = stimOn_times((iChunk)*chunkSize);


unique_templates = unique(spike_templates);
%allUnits_raster = nan(size(allGoodCells,1), size([-0.5:0.1:trialStop-trialStart],2)-1);
allUnits_rasterx = [];
allUnits_rastery = [];
for iCell = 1:size(allGoodCells,1)
    % get spikes 
    raster_window = [-0.5, trialStop-trialStart];
    psth_bin_size = 0.1;
    [~, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(allGoodCells(iCell)), raster_window, psth_bin_size, ...
                        trialStart, []);
    allUnits_rasterx = [allUnits_rasterx , raster_x];
    allUnits_rastery = [allUnits_rastery , raster_y*iCell];
end

stimOn_times_chunk = stimOn_times((1+(iChunk-1)*chunkSize):iChunk*chunkSize)-trialStart;
wheel_move_times_chunk = stim_to_move((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk = stim_to_feedback((1+(iChunk-1)*chunkSize):iChunk*chunkSize).*trial_outcome((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk(feedback_times_chunk==0)=NaN;
trial_cond_chunk = trial_conditions((1+(iChunk-1)*chunkSize):iChunk*chunkSize,:);

figure();

subplot(10,1,1:7)
scatter(t(allUnits_rasterx),  allUnits_rastery, 1, [0,0,0], 'filled')
line([-0.7,-0.7], [1, size(allGoodMSNs,1)], 'Color', [rgb('HotPink'), 1], 'LineWidth', 3);%msn
line([-0.7,-0.7], [size(allGoodMSNs,1)+1, size(allGoodMSNs,1)+size(allGoodTANs,1)],...
    'Color', [rgb('Brown'), 1], 'LineWidth', 3);%tan
line([-0.7,-0.7], [size(allGoodMSNs,1)+size(allGoodTANs,1)+1, ...
    size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1)], 'Color', [rgb('LightSlateGray'), 1], 'LineWidth', 3);%fsi
 
% stimLines = arrayfun(@(x) line([stimOn_times_chunk(x), stimOn_times_chunk(x)],...
%     [1,size(allGoodCells,1)], 'Color', [rgb('Red'), 1], 'LineWidth', 1), 1:size(stimOn_times_chunk,1));
correctGo_times_chunk = stimOn_times_chunk(trial_cond_chunk(:,1)==1 | trial_cond_chunk(:,1)==2 & trial_cond_chunk(:,3)==1);
incorrectGo_times_chunk = stimOn_times_chunk(trial_cond_chunk(:,1)==1 | trial_cond_chunk(:,1)==2 & trial_cond_chunk(:,3)==0);
correctNoGo_times_chunk = stimOn_times_chunk(trial_cond_chunk(:,1)==3 & trial_cond_chunk(:,3)==1);
incorrectNoGo_times_chunk = stimOn_times_chunk(trial_cond_chunk(:,1)==3 & trial_cond_chunk(:,3)==0);

correctGoLines = arrayfun(@(x) line([correctGo_times_chunk(x), correctGo_times_chunk(x)],...
    [1,size(allGoodCells,1)], 'Color', [rgb('DarkGreen'), 1], 'LineWidth', 1), 1:size(correctGo_times_chunk,1));
incorrectGoLines = arrayfun(@(x) line([incorrectGo_times_chunk(x), incorrectGo_times_chunk(x)],...
    [1,size(allGoodCells,1)], 'Color', [rgb('PaleGreen'), 1], 'LineWidth', 1), 1:size(incorrectGo_times_chunk,1));

correctNoGoLines = arrayfun(@(x) line([correctNoGo_times_chunk(x), correctNoGo_times_chunk(x)],...
    [1,size(allGoodCells,1)], 'Color', [rgb('DarkRed'), 1], 'LineWidth', 1), 1:size(correctNoGo_times_chunk,1));
incorrectNoGoLines = arrayfun(@(x) line([incorrectNoGo_times_chunk(x), incorrectNoGo_times_chunk(x)],...
    [1,size(allGoodCells,1)], 'Color', [rgb('Salmon'), 1], 'LineWidth', 1), 1:size(incorrectNoGo_times_chunk,1));

moveLines = arrayfun(@(x) line([stimOn_times_chunk(x) + wheel_move_times_chunk(x), stimOn_times_chunk(x) + wheel_move_times_chunk(x)],...
    [size(allGoodCells,1)+5, size(allGoodCells,1)], 'Color', [rgb('Purple'), 1], 'LineWidth', 1), 1:size(wheel_move_times_chunk,1));
rewLines = arrayfun(@(x) line([stimOn_times_chunk(x) + feedback_times_chunk(x), stimOn_times_chunk(x) + feedback_times_chunk(x)],...
    [size(allGoodCells,1)+5, size(allGoodCells,1)], 'Color', [rgb('Blue'), 1], 'LineWidth', 1), 1:size(feedback_times_chunk ,1));
xlim([-0.7, t(end)])
xlabel('time (s)')
ylabel('cell #')

subplot(10,1,8)
plot(trial_conditions((1+(iChunk-1)*chunkSize):iChunk*chunkSize,1)); 
ylabel('stim type')

subplot(10,1,8)%correct go 
plot(trial_conditions((1+(iChunk-1)*chunkSize):iChunk*chunkSize,1)); 

subplot(10,1,8)%correct no go 
plot(wheel_velocity((1+(iChunk-1)*chunkSize):iChunk*chunkSize,1)); 


% psth all incorrect, all correct ect... 
%cell * time, avarge accross trials 

%allUnits_raster = nan(size(allGoodCells,1), size([-0.5:0.1:trialStop-trialStart],2)-1);
correctGo_early =  (trial_conditions(:,1)==2 | trial_conditions(:,1)==1) & trial_conditions(:,3)==1 &stim_to_move <=0.22;
correctGo_late =  (trial_conditions(:,1)==2 | trial_conditions(:,1)==1) & trial_conditions(:,3)==1 &stim_to_move >0.22;
incorrectGo = (trial_conditions(:,1)==2 | trial_conditions(:,1)==1) & trial_conditions(:,3)==0;

correctNoGo = (trial_conditions(:,1)==3) & trial_conditions(:,3)==1;
incorrectNoGo_early =  (trial_conditions(:,1)==3) & trial_conditions(:,3)==0&stim_to_move <=0.22;
incorrectNoGo_late =  (trial_conditions(:,1)==3) & trial_conditions(:,3)==0&stim_to_move >0.22;

trialType = nan(size(correctGo,1),1);
trialType(correctGo_early) = 1;
trialType(correctGo_late) = 2;
trialType(incorrectGo) = 3;
trialType(correctNoGo) = 4;
trialType(incorrectNoGo_early) = 5;
trialType(incorrectNoGo_late) = 6;

cell_psth=nan(size(allGoodCells,1),100,6);
cell_psth_all=nan(size(allGoodCells,1),100);

for iCell = 1:size(allGoodCells,1)
    raster_window = [-0.5, 0.5];
    psth_bin_size = 0.01;
    [curr_psth, ~, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(allGoodCells(iCell)), raster_window, psth_bin_size, ...
                        stimOn_times, trialType);
    cell_psth(iCell,:,:) = curr_psth';

     [curr_psth_all, ~, t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(allGoodCells(iCell)), raster_window, psth_bin_size, ...
                        stimOn_times, []);
    cell_psth_all(iCell,:) = curr_psth_all;
end

% zscore cell psth 
zscore_correct_go_early = (cell_psth(:,:,1) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
figure();
subplot(161)
imagesc(t,[],zscore_correct_go_early); hold on;
line([-0.5,-0.5], [0.5, size(allGoodMSNs,1)+0.5], 'Color', [rgb('HotPink'), 1], 'LineWidth', 3);%msn
line([-0.5,-0.5], [size(allGoodMSNs,1)+0.5, size(allGoodMSNs,1)+size(allGoodTANs,1)+0.5],...
    'Color', [rgb('Brown'), 1], 'LineWidth', 3);%tan
line([-0.5,-0.5], [size(allGoodMSNs,1)+size(allGoodTANs,1)+0.5, ...
    size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1)+0.5], 'Color', [rgb('LightSlateGray'), 1], 'LineWidth', 3);%fsi
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
xlim([-0.52, 0.5])
%colorbar
title('Correct go early')

subplot(162)
zscore_correct_go_late = (cell_psth(:,:,2) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
imagesc(t,[],zscore_correct_go_late); hold on;
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
%xlim([-0.7, 0.5])
%colorbar
title('Correct go late')

zscore_incorrect_go = (cell_psth(:,:,3) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
subplot(163)
imagesc(t,[],zscore_incorrect_go)
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
%colorbar
title('Incorrect nogo')

zscore_correct_nogo = (cell_psth(:,:,4) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
subplot(164)
imagesc(t,[],zscore_correct_nogo)
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
%colorbar
title('Correct nogo')

zscore_incorrect_nogo_early = (cell_psth(:,:,5) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
subplot(165)
imagesc(t,[],zscore_incorrect_nogo)
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
%colorbar
title('Incorrect go early')

zscore_incorrect_nogo_late = (cell_psth(:,:,6) - nanmean(cell_psth_all(:,1:50),2)) ./ nanstd(cell_psth_all(:,1:50),[],2);
subplot(166)
imagesc(t,[],zscore_incorrect_nogo_late)
colormap(brewermap([],'*RdBu'));
clim([-20, 20])
colorbar
title('Incorrect go late')

% distribution recation times? 
figure();
subplot(311)
histogram(stim_to_move(correctGo), 0:0.02:1.8)
ylabel('count #')
title('reaction time, correct go')
ylims=ylim;
line([0.22, 0.22], [0, ylims(2)], 'Color', 'red')

subplot(312)
histogram(stim_to_move(incorrectNoGo), 0:0.02:1.8)
ylabel('count #')
title('reaction time, incorrect go')
ylims=ylim;
line([0.22, 0.22], [0, ylims(2)], 'Color', 'red')

subplot(313)
histogram(stim_to_move(incorrectGo), 0:0.02:1.8)
ylabel('count #')
xlabel('time from stim onset(s)')
title('reaction time, incorrect nogo')
ylims=ylim;
line([0.22, 0.22], [0, ylims(2)], 'Color', 'red')

%% scatter average baseline f.r. two plots correct go/incorrect go ; and for nogo 
figure();
subplot(121)
scatter(nanmean(cell_psth(1:size(allGoodMSNs,1),1:50,1),2).*50, nanmean(cell_psth(1:size(allGoodMSNs,1),1:50,3),2).*50, 9, rgb('HotPink'),'filled'); hold on;
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),1:50,1),2).*50,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),1:50,3),2).*50, 15, rgb('Brown'),'filled')
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),1:50,1),2).*50,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),1:50,3),2).*50, 12, rgb('SlateGrey'),'filled')
xlabel('f.r. correct go')
ylabel('f.r. incorrect nogo')
xlim([0.01, 14])
ylim([0.01, 14])
line([0.01, 14], [0.01, 14], 'Color', 'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

subplot(122)
scatter(nanmean(cell_psth(1:size(allGoodMSNs,1),1:50,4),2).*50, nanmean(cell_psth(1:size(allGoodMSNs,1),1:50,5),2).*50, 9, rgb('HotPink'),'filled'); hold on;
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),1:50,4),2).*50,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),1:50,5),2).*50, 15, rgb('Brown'),'filled')
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),1:50,4),2).*50,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),1:50,5),2).*50, 12, rgb('SlateGrey'),'filled')
xlabel('f.r. correct nogo')
ylabel('f.r. incorrect go')
xlim([0.01, 14])
ylim([0.01, 14])
line([0.01, 14], [0.01, 14], 'Color', 'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


figure();
subplot(121)
scatter(nanmean(cell_psth(1:size(allGoodMSNs,1),55:69,2),2).*14, nanmean(cell_psth(1:size(allGoodMSNs,1),55:69,3),2).*14, 9, rgb('HotPink'),'filled'); hold on;
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),55:69,2),2).*14,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),55:69,3),2).*14, 15, rgb('Brown'),'filled')
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),55:69,2),2).*14,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),55:69,3),2).*14, 12, rgb('SlateGrey'),'filled')
xlabel('f.r. correct go')
ylabel('f.r. incorrect nogo')
xlim([0.01, 14])
ylim([0.01, 14])
line([0.01, 14], [0.01, 14], 'Color', 'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

subplot(122)
scatter(nanmean(cell_psth(1:size(allGoodMSNs,1),55:69,4),2).*14, nanmean(cell_psth(1:size(allGoodMSNs,1),55:69,5),2).*14, 9, rgb('HotPink'),'filled'); hold on;
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),55:69,4),2).*14,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1:size(allGoodMSNs,1)+size(allGoodTANs,1),55:69,5),2).*14, 15, rgb('Brown'),'filled')
scatter(nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),55:69,4),2).*14,...
    nanmean(cell_psth(size(allGoodMSNs,1)+1+size(allGoodTANs,1):size(allGoodMSNs,1)+size(allGoodTANs,1)+size(allGoodFSIs,1),55:69,5),2).*14, 12, rgb('SlateGrey'),'filled')
xlabel('f.r. correct nogo')
ylabel('f.r. incorrect go')
xlim([0.01, 14])
ylim([0.01, 14])
line([0.01, 14], [0.01, 14], 'Color', 'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% correlation pupil size lapse rate 
% scatter average stim response 