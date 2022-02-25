
function psthGUI(spike_templates, task_spike_timeline, task_stimOn_times, task_wheel_starts, ...
     task_reward,  task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
    passive_stimOn_times, passive_wheel_starts )

%% set up dynamic figure
psthGuiHandle = figure('color', 'w');
set(psthGuiHandle, 'KeyPressFcn', @KeyPressCb);

%% initial conditions
iCluster = 1;
uniqueTemps = unique(spike_templates);

%% plot initial conditions
initializePlot(psthGuiHandle)
updateUnit(iCluster, uniqueTemps,  spike_templates, psthGuiHandle,task_spike_timeline, task_stimOn_times, task_wheel_starts, ...
    task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
    passive_stimOn_times, passive_wheel_starts);

%% change on keypress
    function KeyPressCb(~, evnt)
        %fprintf('key pressed: %s\n', evnt.Key);
        if strcmpi(evnt.Key, 'rightarrow')
            iCluster = iCluster + 1;
            updateUnit(iCluster, uniqueTemps,  spike_templates, psthGuiHandle,task_spike_timeline, task_stimOn_times,  task_wheel_starts, ...
                 task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
                passive_stimOn_times, passive_wheel_starts);
        elseif strcmpi(evnt.Key, 'leftarrow')
            iCluster = iCluster - 1;
            updateUnit(iCluster, uniqueTemps,  spike_templates, psthGuiHandle,task_spike_timeline, task_stimOn_times,  task_wheel_starts, ...
                 task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
                passive_stimOn_times, passive_wheel_starts);
        end
    end
end

function initializePlot(psthGuiHandle)

%% main title

mainTitle = sgtitle('');

%% initialize task stim raster

ss=subplot(6, 2, [1, 3]);
title('task stim')
hold on;
rasterTaskDots = scatter(NaN,NaN,4,'k','filled');
rasterTaskStimLine = line([NaN, NaN], [NaN, NaN], 'Color', 'r');
rasterTaskMoveDots = scatter(NaN,NaN,4,'g','filled');
rasterTaskRewaLine = scatter(NaN,NaN,4,'b','filled');
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')
makepretty;
%% initialize task stim PSTH 
subplot(6, 2, 5);
hold on;
psthTaskLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')
makepretty;
%% initialize passive stim raster
subplot(6, 2, [2, 4])
title('passive stim')
hold on;
rasterPassiveDots = scatter(NaN,NaN,4,'k','filled');
rasterPassiveStimLine = line([NaN, NaN], [NaN, NaN],'Color', 'r');
rasterPassiveMoveDots = scatter(NaN,NaN,4,'g','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')
makepretty;
%% initialize passive stim psth 
subplot(6, 2, 6)
hold on;
psthPassiveLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from stim (s)')
makepretty;
%% initialize task move raster

subplot(6, 2, [7, 9]);
title('task move')
hold on;
rasterTaskMoveMoveDots = scatter(NaN,NaN,4,'k','filled');
rasterTaskMoveLine = line([NaN, NaN], [NaN, NaN],'Color', 'g');
rasterTaskStimMoveDots = scatter(NaN,NaN,4,'r','filled');
rasterTaskRewaMoveLine = scatter(NaN,NaN,4,'b','filled');
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from move (s)')
makepretty;
%% initialize task stim PSTH 
subplot(6, 2, 11);
hold on;
psthTaskMoveLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from move (s)')
makepretty;
%% initialize task move sponaneous stim raster
subplot(6, 2, [8, 10]);
title('task spont move')
hold on;
rasterTaskMoveSpontDots = scatter(NaN,NaN,4,'k','filled');
rasterTaskMoveSpontLine = line([NaN, NaN], [NaN, NaN],'Color', 'g');
rasterTaskStimMoveSpontDots = scatter(NaN,NaN,4,'r','filled');
rasterTaskRewaMoveSpontLine = scatter(NaN,NaN,4,'b','filled');
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from move (s)')
makepretty;
%% initialize task stim PSTH 
subplot(6, 2, 12);
hold on;
psthTaskMoveSpontLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from move (s)')
makepretty;
%% save all handles
guiData = struct;
% main title
guiData.mainTitle = mainTitle;
% task stim raster
guiData.ss=ss;
guiData.rasterTaskDots = rasterTaskDots;
guiData.rasterTaskStimLine = rasterTaskStimLine;
guiData.rasterTaskMoveDots = rasterTaskMoveDots;
guiData.rasterTaskRewaLine = rasterTaskRewaLine;
% task stim psth 
guiData.psthTaskLines = psthTaskLines; 
% passive stim raster
guiData.rasterPassiveDots = rasterPassiveDots;
guiData.rasterPassiveStimLine = rasterPassiveStimLine;
guiData.rasterPassiveMoveDots = rasterPassiveMoveDots;
% passive stim psth 
guiData.psthPassiveLines = psthPassiveLines;
% task move raster
guiData.rasterTaskMoveMoveDots = rasterTaskMoveMoveDots;
guiData.rasterTaskMoveLine = rasterTaskMoveLine;
guiData.rasterTaskStimMoveDots = rasterTaskStimMoveDots;
guiData.rasterTaskRewaMoveLine = rasterTaskRewaMoveLine;
% task move psth
guiData.psthTaskMoveLines = psthTaskMoveLines;
% task move spont raster
guiData.rasterTaskMoveSpontDots = rasterTaskMoveSpontDots;
guiData.rasterTaskMoveSpontLine = rasterTaskMoveSpontLine;
guiData.rasterTaskStimMoveSpontDots = rasterTaskStimMoveSpontDots;
guiData.rasterTaskRewaMoveSpontLine = rasterTaskRewaMoveSpontLine;
% task move spont psth
guiData.psthTaskMoveSpontLines = psthTaskMoveSpontLines;
% upload guiData
guidata(psthGuiHandle, guiData);
end

function updateUnit(iCluster, uniqueTemps,  spike_templates, psthGuiHandle,task_spike_timeline, task_stimOn_times,  task_wheel_starts, ...
                 task_reward, task_stim_to_move, task_stim_to_feedback, passive_spike_timeline, ...
                passive_stimOn_times, passive_wheel_starts)

%% Get guidata
guiData = guidata(psthGuiHandle);
thisUnit = uniqueTemps(iCluster);
set(guiData.mainTitle, 'String', num2str(thisUnit))
[raster_x, raster_y, t, curr_smoothed_psth,trial_sort]= psthGet(spike_templates,thisUnit,task_spike_timeline,task_stimOn_times, task_stim_to_move);
%% plot task stim raster

set(guiData.rasterTaskDots,'XData',t(raster_x),'YData',raster_y);
set(guiData.rasterTaskStimLine, 'XData', [0, 0], 'YData', [min(raster_y), max(raster_y)])
set(guiData.rasterTaskMoveDots, 'XData', task_stim_to_move(trial_sort), 'YData', 1:numel(trial_sort))
set(guiData.rasterTaskRewaLine, 'XData', task_stim_to_feedback(trial_sort), 'YData', 1:numel(trial_sort))
set(guiData.ss, 'xlim',[-0.3, 0.5]);
%% plot task stim psth
set(guiData.psthTaskLines, ...
    'XData',t,'YData',curr_smoothed_psth)
%% plot passive stim raster 
[raster_x, raster_y, t, curr_smoothed_psth,trial_sort,curr_raster_sorted]= psthGet(spike_templates,thisUnit,passive_spike_timeline,passive_stimOn_times, ones(length(passive_stimOn_times),1));

set(guiData.rasterPassiveDots,'XData',t(raster_x),'YData',raster_y);
set(guiData.rasterPassiveStimLine, 'XData', [0, 0], 'YData', [min(raster_y), max(raster_y)])
%% plot passive stim psth
frB = sum(curr_raster_sorted(1:2:end,100:290));
frA = sum(curr_raster_sorted(1:2:end,291:481));
p=signrank(frB,frA);
frB = sum(curr_raster_sorted(2:2:end,100:290));
frA = sum(curr_raster_sorted(2:2:end,291:481));
p2=signrank(frB,frA);
if p<0.05 && p2 <0.05
    set(guiData.psthPassiveLines, ...
    'XData',t,'YData',curr_smoothed_psth, 'CData', rgb('Orange'))
else
    set(guiData.psthPassiveLines, ...
    'XData',t,'YData',curr_smoothed_psth,'CData', rgb('Black'))
end

[raster_x, raster_y, t, curr_smoothed_psth, trial_sort,curr_raster_sorted]= psthGet(spike_templates,thisUnit,task_spike_timeline,task_wheel_starts, ones(length(task_wheel_starts),1));
%% plot task move raster
set(guiData.rasterTaskMoveSpontDots,'XData',t(raster_x),'YData',raster_y);
set(guiData.rasterTaskMoveSpontLine, 'XData', [0, 0], 'YData', [min(raster_y), max(raster_y)])

%% plot task move psth

frB = sum(curr_raster_sorted(1:2:end,100:290));
frA = sum(curr_raster_sorted(1:2:end,391:581));
p=signrank(frB,frA);
frB = sum(curr_raster_sorted(2:2:end,100:290));
frA = sum(curr_raster_sorted(2:2:end,391:581));
p2=signrank(frB,frA);
if p<0.05 && p2 <0.05
    set(guiData.psthTaskMoveSpontLines, ...
    'XData',t,'YData',curr_smoothed_psth, 'CData', rgb('Orange'))
else
    set(guiData.psthTaskMoveSpontLines, ...
    'XData',t,'YData',curr_smoothed_psth,'CData', rgb('Black'))
end
[raster_x, raster_y, t, curr_smoothed_psth]= psthGet(spike_templates,thisUnit,task_spike_timeline,task_stimOn_times + task_stim_to_move, ones(length(task_stimOn_times),1));
%% plot task move raster
set(guiData.rasterTaskMoveMoveDots,'XData',t(raster_x),'YData',raster_y);
set(guiData.rasterTaskMoveLine, 'XData', [0, 0], 'YData', [min(raster_y), max(raster_y)])

%% plot task move psth
set(guiData.psthTaskMoveLines, ...
    'XData',t,'YData',curr_smoothed_psth)

end


function [raster_x, raster_y, t, curr_smoothed_psth,trial_sort,curr_raster_sorted]= psthGet(spike_templates,thisUnit,spike_timeline, thisAlign, thisSort)
raster_window = [-0.3,0.5];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
t_peri_event = thisAlign + t_bins;
% (handle NaNs by setting rows with NaN times to 0)
t_peri_event(any(isnan(t_peri_event),2),:) = 0;

% Bin spikes (use only spikes within time range, big speed-up)
curr_spikes_idx = ismember(spike_templates,thisUnit);
curr_raster_spike_times = spike_timeline(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
    curr_raster_spike_times > max(t_peri_event(:))) = [];

if ~any(diff(reshape(t_peri_event',[],1)) < 0)
    % (if no backward time jumps, can do long bin and cut out in-between, faster)
    curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
        reshape(t_peri_event',[],1)),NaN],size(t_peri_event'))';
    curr_raster = curr_raster_continuous(:,1:end-1);   
else
    % (otherwise, bin trial-by-trial)
    curr_raster = cell2mat(arrayfun(@(x) ...
        histcounts(curr_raster_spike_times,t_peri_event(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false));
end

smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_t = mean(diff(t));

curr_psth =  mean(curr_raster,1);
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;
[~,trial_sort] = sort(thisSort);
curr_raster_sorted = curr_raster(trial_sort,:);
[raster_y,raster_x] = find(curr_raster_sorted);
end