
function psthGUI_precalc(rasterX_stim, rasterY_stim, PSTHstim, rasterXmove, rasterYmove, PSTHmove, t, t_move)

%% set up dynamic figure
psthGuiHandle = figure('color', 'w');
set(psthGuiHandle, 'KeyPressFcn', @KeyPressCb);

%% initial conditions
iCluster = 1;

%% plot initial conditions
initializePlot(psthGuiHandle)
updateUnit(iCluster, psthGuiHandle, rasterX_stim, rasterY_stim, PSTHstim, rasterXmove, rasterYmove, PSTHmove,t,t_move);

%% change on keypress
    function KeyPressCb(~, evnt)
        %fprintf('key pressed: %s\n', evnt.Key);
        if strcmpi(evnt.Key, 'rightarrow')
            iCluster = iCluster + 1;
            updateUnit(iCluster, psthGuiHandle, rasterX_stim, rasterY_stim, PSTHstim, rasterXmove, rasterYmove, PSTHmove,t,t_move);
        elseif strcmpi(evnt.Key, 'leftarrow')
            iCluster = iCluster - 1;
            updateUnit(iCluster, psthGuiHandle, rasterX_stim, rasterY_stim, PSTHstim, rasterXmove, rasterYmove, PSTHmove,t,t_move);
        end
    end
end

function initializePlot(psthGuiHandle)

%% main title

mainTitle = sgtitle('');

%% initialize task stim raster
% 
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
ylabel('Sp/s')
xlabel('time from stim (s)')
makepretty;
%% initialize passive stim raster
subplot(6, 2, [2, 4])
%title('passive stim')
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
ylabel('Sp/s')
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
%% initialize task move PSTH 
subplot(6, 2, 11);
hold on;
psthTaskMoveLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('Sp/s')
xlabel('time from move (s)')
makepretty;
%% initialize task move sponaneous stim raster
subplot(6, 2, [8, 10]);
%title('task spont move')
hold on;
rasterTaskMoveSpontDots = scatter(NaN,NaN,4,'k','filled');
rasterTaskMoveSpontLine = line([NaN, NaN], [NaN, NaN],'Color', 'g');
rasterTaskStimMoveSpontDots = scatter(NaN,NaN,4,'r','filled');
rasterTaskRewaMoveSpontLine = scatter(NaN,NaN,4,'b','filled');
xlim([-0.3, 0.5]);
ylabel('trial #')
xlabel('time from move (s)')
makepretty;
%% initialize task move sponaneous stim PSTH 
subplot(6, 2, 12);
hold on;
psthTaskMoveSpontLines = scatter(NaN,NaN,4,'k','filled');
%psthTaskStimLine = 
xlim([-0.3, 0.5]);
ylabel('Sp/s')
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

function updateUnit(iCluster, psthGuiHandle, rasterX_stim, rasterY_stim, PSTHstim, rasterX_move, rasterY_move, PSTHmove,t, t_move)

%% Get guidata
guiData = guidata(psthGuiHandle);
thisUnit = iCluster;
set(guiData.mainTitle, 'String', num2str(thisUnit))
%% plot task stim raster
% set(guiData.rasterTaskDots,'XData',t(rasterX_stim),'YData',rasterY_stim);
% set(guiData.rasterTaskStimLine, 'XData', [0, 0], 'YData', [min(rasterY_stim), max(rasterY_stim)])
% set(guiData.rasterTaskMoveDots, 'XData', task_stim_to_move(trial_sort), 'YData', 1:numel(trial_sort))
% set(guiData.rasterTaskRewaLine, 'XData', task_stim_to_feedback(trial_sort), 'YData', 1:numel(trial_sort))
% set(guiData.ss, 'xlim',[-0.3, 0.5]);
% %% plot task stim psth
% set(guiData.psthTaskLines, ...
%     'XData',t,'YData',PSTHstim)
%% plot passive stim raster 

set(guiData.rasterPassiveDots,'XData',t(rasterX_stim{thisUnit}),'YData',rasterY_stim{thisUnit});
set(guiData.rasterPassiveStimLine, 'XData', [0, 0], 'YData', [min(rasterY_stim{thisUnit}), max(rasterY_stim{thisUnit})])
set(guiData.ss, 'xlim',[-0.3, 0.5]);
%% plot passive stim psth

set(guiData.psthPassiveLines, ...
    'XData',t,'YData',PSTHstim(thisUnit,:), 'CData', rgb('Black'))


%% plot task move raster
set(guiData.rasterTaskMoveSpontDots,'XData',t_move(rasterX_move{thisUnit}),'YData',rasterY_move{thisUnit});
set(guiData.rasterTaskMoveSpontLine, 'XData', [0, 0], 'YData', [min(rasterY_move{thisUnit}), max(rasterY_move{thisUnit})])

%% plot task move psth


    set(guiData.psthTaskMoveSpontLines, ...
    'XData',t_move,'YData',PSTHmove(thisUnit,:), 'CData', rgb('Black'))

%% plot task move raster
% set(guiData.rasterTaskMoveMoveDots,'XData',t(raster_x),'YData',raster_y);
% set(guiData.rasterTaskMoveLine, 'XData', [0, 0], 'YData', [min(raster_y), max(raster_y)])

%% plot task move psth
% set(guiData.psthTaskMoveLines, ...
%     'XData',t_move,'YData',PSTHmove(thisUnit))

end




