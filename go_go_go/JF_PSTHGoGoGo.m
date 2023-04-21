
%% population PSTH
% load all data in one region, get psth per cell (baseline z-scored)

%urr_psth_stim_active = nan(1,4,100);
clearvars curr_psth_stim_active curr_psth_rew_active curr_psth_stim_passive curr_psth_move_active
unitCount = 0;
for iRecording = 1:size(ephysData, 2)
    uniqueTemps = ephysData(iRecording).spike_templates;

    for iUnit = 1:size(ephysData(iRecording).unitType, 1)
        unitCount = unitCount+1;
        thisUnit = uniqueTemps(iUnit);
        unitTypes(unitCount) = ephysData(iRecording).unitType(iUnit);
        prop_isi(unitCount) = ephysData(iRecording).prop_isi(iUnit);
        wav_dur(unitCount) = ephysData(iRecording).waveform_duration(iUnit);
        % aligned to stim onset per image in active (late move trials)
        late_move_trails = isnan(ephysData(iRecording).stim_to_move) | ephysData(iRecording).stim_to_move >= 0.35;
        threeStims = ismember(ephysData(iRecording).trial_conditions(:,1), 1:3);
        [~, curr_psth_stim_active(unitCount, :, :), ~, ~, ~, ~] = JF_raster_PSTH(ephysData(iRecording).spike_templates, ...
            ephysData(iRecording).spike_times, thisUnit, [-0.5, 0.35], 0.001, ...
            ephysData(iRecording).stimOn_times(late_move_trails & threeStims), ...
            ephysData(iRecording).trial_conditions(late_move_trails & threeStims, 1), ...
            [], [], 0, 1, []);

        % aligned to stim onset per image in passive (no move trials)
        theseImages = [4,6,12];
        late_move_trails = isnan(ephysData(iRecording).stim_to_move_passive) | ephysData(iRecording).stim_to_move_passive >= 0.35;
        threeStims = ismember(ephysData(iRecording).trial_conditions_passive(:,1), theseImages);
        controlateral = ephysData(iRecording).trial_conditions_passive(:,2) == -90;
        [~, curr_psth_stim_passive(unitCount, :, :), ~, ~, ~, ~] = JF_raster_PSTH(ephysData(iRecording).spike_templates, ...
            ephysData(iRecording).spike_times, thisUnit, [-0.5, 0.35], 0.001, ...
            ephysData(iRecording).stimOn_times_passive(late_move_trails & threeStims & controlateral), ...
            ephysData(iRecording).trial_conditions_passive(late_move_trails & threeStims & controlateral, 1), ...
            [], [], 0, 1, []);

        % aligned to move onset per move
        [~, curr_psth_move_active(unitCount, :, :), ~, ~, ~, ~] = JF_raster_PSTH(ephysData(iRecording).spike_templates, ...
            ephysData(iRecording).spike_times, thisUnit, [-0.5, 0.5], 0.001, ...
            ephysData(iRecording).stimOn_times + ephysData(iRecording).stim_to_move, ...
            ones(size(ephysData(iRecording).trial_conditions(:,2),1),1), ...
            [], [], 0, 1, []); %QQ not distinguishing movements for now b/c some recordings only have one type (and very few wrong ones in general) 

        % aligned to reward onset per reward type
        [~, curr_psth_rew_active(unitCount, :, :), ~, ~, ~, ~] = JF_raster_PSTH(ephysData(iRecording).spike_templates, ...
            ephysData(iRecording).spike_times, thisUnit, [-0.5, 0.5], 0.001, ...
            ephysData(iRecording).rewards_times, ...
            ephysData(iRecording).rewards_type, ...
            [], [], 0, 1, []);
    end


end

% classify cells
msns = prop_isi < 30 & wav_dur > 400;
fsi = prop_isi < 30 & wav_dur <= 400;
tan = prop_isi >= 30;
cell_types = [msns + (fsi*2) + (tan *3)];

% get stim response 
figure();
for iCellType =1:3
    stim_psth = squeeze(curr_psth_stim_passive(unitTypes==1 & cell_types == iCellType, 1, :));
    %stim_psth_zscore = (stim_psth - nanmean(stim_psth(:,1:500),2)) ./  nanstd(stim_psth(:,1:500),[],2);
    smooth_filt = [ceil(size(stim_psth,1)./40),20]; % (units x frames)
    
    
    %max−min/value−min
    stim_psth_max_min = (stim_psth- min(stim_psth,[],2))./...
        (min(stim_psth,[],2)+max(stim_psth,[],2));
     [unitSort, unitSort_idx] = sort(nanmean(stim_psth_max_min(:,550:750),2));
     stim_psth_zscore_smooth = conv2(stim_psth_max_min(unitSort_idx,:),ones(smooth_filt),'same')./ ...
                conv2(~isnan(stim_psth_max_min(unitSort_idx,:)),ones(smooth_filt),'same');
   
    subplot(1,3,iCellType)
    imagesc(-0.5:0.01:0.34, [],stim_psth_zscore_smooth)
    colormap(brewermap([],'*RdBu'));
    clim([-abs(max(max(stim_psth_zscore_smooth))), abs(max(max(stim_psth_zscore_smooth)))])
    colorbar;
end

% global psths
% zscore the traces
rew_psth = squeeze(curr_psth_rew_active(unitTypes==1 & prop_isi > 30, 2, :));
rew_psth_zscore = (rew_psth - nanmean(rew_psth(:,1:500),2)) ./  nanstd(rew_psth(:,1:500),[],2);
[unitSort, unitSort_idx] = sort(nanmean(rew_psth_zscore(:,550:750),2));
smooth_filt = [10,15]; % (units x frames)

rew_psth_zscore_smooth = conv2(rew_psth_zscore(unitSort_idx,:),ones(smooth_filt),'same')./ ...
            conv2(~isnan(rew_psth_zscore(unitSort_idx,:)),ones(smooth_filt),'same');

figure();
imagesc(-0.5:0.01:0.34, [],rew_psth_zscore_smooth)
colormap(brewermap([],'*RdBu'));
clim([-abs(max(max(rew_psth_zscore_smooth))), abs(max(max(rew_psth_zscore_smooth)))])
colorbar;



stim_psth = squeeze(curr_psth_stim_active(unitTypes==1, 1, :));
stim_psth_zscore = (stim_psth - nanmean(stim_psth(:,1:500),2)) ./  nanstd(stim_psth(:,1:500),[],2);
[unitSort, unitSort_idx] = sort(nanmean(stim_psth_zscore(:,550:750),2));
smooth_filt = [10,15]; % (units x frames)

stim_psth_zscore_smooth = conv2(stim_psth_zscore(unitSort_idx,:),ones(smooth_filt),'same')./ ...
            conv2(~isnan(stim_psth_zscore(unitSort_idx,:)),ones(smooth_filt),'same');

figure();
imagesc(-0.5:0.01:0.34, [],stim_psth_zscore_smooth)
colormap(brewermap([],'*RdBu'));
clim([-abs(max(max(stim_psth_zscore_smooth))), abs(max(max(stim_psth_zscore_smooth)))])
colorbar;

figure();
stim_psth = squeeze(curr_psth_stim_active(unitTypes==1, 2, :));
stim_psth_zscore = (stim_psth - nanmean(stim_psth(:,1:500),2)) ./  nanstd(stim_psth(:,1:500),[],2);
smooth_filt = [10,15]; % (units x frames)

stim_psth_zscore_smooth = conv2(stim_psth_zscore(unitSort_idx,:),ones(smooth_filt),'same')./ ...
            conv2(~isnan(stim_psth_zscore(unitSort_idx,:)),ones(smooth_filt),'same');

imagesc(-0.5:0.01:0.34, [],stim_psth_zscore_smooth)
colormap(brewermap([],'*RdBu'));
clim([-abs(max(max(stim_psth_zscore_smooth))), abs(max(max(stim_psth_zscore_smooth)))])
colorbar;
% cell classification

%% summary

%% trial PSTH

%% example cells