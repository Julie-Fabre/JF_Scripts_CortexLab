iChunk = 1;
chunkSize = 454;
trialStart = stimOn_times(1+(iChunk-1)*chunkSize);
trialStop = stimOn_times((iChunk)*chunkSize);


unique_templates = unique(spike_templates);
allUnits_raster = nan(size(allGoodCells,1), size([-0.5:0.1:trialStop-trialStart],2)-1);
allUnits_rasterx = [];
allUnits_rastery = [];
for iCell = 1:size(allGoodCells,1)
    % get spikes 
    raster_window = [-0.5, trialStop-trialStart];
    psth_bin_size = 0.1;
    [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(allGoodCells(iCell)), raster_window, psth_bin_size, ...
                        trialStart, []);
    allUnits_rasterx = [allUnits_rasterx , raster_x];
    allUnits_rastery = [allUnits_rastery , raster_y*iCell];
    allUnits_raster(iCell,:) = curr_raster;
end

stimOn_times_chunk = stimOn_times((1+(iChunk-1)*chunkSize):iChunk*chunkSize)-trialStart;
wheel_move_times_chunk = stim_to_move((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk = stim_to_feedback((1+(iChunk-1)*chunkSize):iChunk*chunkSize).*trial_outcome((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk(feedback_times_chunk==0)=NaN;
trial_cond_chunk = trial_conditions((1+(iChunk-1)*chunkSize):iChunk*chunkSize,:);

figure();hold on;
plot(allUnits_raster(:,1:1000)')


%% rastermap
ops.nCall = [5, 50]
[isort1, isort2, Sm] = mapTmap(allUnits_raster, ops);
imagesc(allUnits_raster(isort1,:))

%The matlab code needs to be cleaned up but the main function to call is mapTmap.m. This function is used in the example script loadFromPython.m (loads suite2p outputs, requires npy-matlab).
[iclustup, isort, Vout] = activityMap(allUnits_raster);
figure();
imagesc(allUnits_raster(isort,:))

% ok that doens't seem to really be working... let's just sort by firing
% rate to stim 
[~, frGo_idx_msn] = sort((nanmean(cell_psth(1:size(allGoodMSNs,1),55:69,2),2)-nanmean(cell_psth(1:size(allGoodMSNs,1),1:50,2),2))./nanstd(cell_psth(1:size(allGoodMSNs,1),1:50,2),[],2));
[~, frGo_idx_tan] = sort((nanmean(cell_psth(1:size(allGoodTANs,1),55:69,2),2)-nanmean(cell_psth(1:size(allGoodTANs,1),1:50,2),2))./nanstd(cell_psth(1:size(allGoodTANs,1),1:50,2),[],2));
[~, frGo_idx_fsi] = sort((nanmean(cell_psth(1:size(allGoodFSIs,1),55:69,2),2)-nanmean(cell_psth(1:size(allGoodFSIs,1),1:50,2),2))./nanstd(cell_psth(1:size(allGoodFSIs,1),1:50,2),[],2));
fr_idx = [allGoodCells(frGo_idx_msn); allGoodCells(size(allGoodMSNs,1)+frGo_idx_tan); allGoodCells(size(allGoodMSNs,1)+size(allGoodTANs,1)+frGo_idx_fsi)];
%% trials
iChunk = iChunk+1;
chunkSize = 10;
trialStart = stimOn_times(1+(iChunk-1)*chunkSize);
trialStop = stimOn_times((iChunk)*chunkSize);


unique_templates = unique(spike_templates);
allUnits_raster = nan(size(allGoodCells,1), size([-0.5:0.1:trialStop-trialStart],2)-1);
allUnits_rasterx = [];
allUnits_rastery = [];
for iCell = 1:size(allGoodCells,1)
    % get spikes 
    raster_window = [-0.5, trialStop-trialStart];
    psth_bin_size = 0.001;
    [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
                        unique_templates(fr_idx(iCell)), raster_window, psth_bin_size, ...
                        trialStart, []);
    allUnits_rasterx = [allUnits_rasterx , raster_x];
    allUnits_rastery = [allUnits_rastery , raster_y*iCell];
    allUnits_raster(iCell,:) = curr_psth;
end

stimOn_times_chunk = stimOn_times((1+(iChunk-1)*chunkSize):iChunk*chunkSize)-trialStart;
wheel_move_times_chunk = stim_to_move((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk = stim_to_feedback((1+(iChunk-1)*chunkSize):iChunk*chunkSize).*trial_outcome((1+(iChunk-1)*chunkSize):iChunk*chunkSize);
feedback_times_chunk(feedback_times_chunk==0)=NaN;
trial_cond_chunk = trial_conditions((1+(iChunk-1)*chunkSize):iChunk*chunkSize,:);

clf;
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