
%% sanity check task 
clear all;

animal = 'JF093';
day = '2023-03-06';
% go no go example: animal JF093, day 2023-03-06 exp 2
% go go go example: JF096, day 2023-02-07 exp 2
experiment = 1;
load_parts.cam = true;
verbose = true;
JF_load_experiment;
camera_location = facecam_fn;
camera_t = facecam_t;
% Get wheel movements during stim, only use quiescent trials

framerate = 30;

wheel_window = [0, 0.5];
wheel_window_t = wheel_window(1):1 / framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus, stimOn_times, wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0, 2);

timeline_rate = 0.0010;

stims = [1,2,3];
stim_types = {'Stim 1', 'Stim 2', 'Stim 3'};
%surround_frames1 = 30;
surround_frames = 150;
surround_t = [-surround_frames:surround_frames]./framerate;


%use_stim = stims(iStim);
use_align = stimOn_times;% + stim_to_move(~isnan(stim_to_move));%stimOn_times(stimIDs == use_stim & quiescent_trials & no_move_trials==1);

% Initialize video reader, get average and average difference
box_crop = [];
plot_me = 0;
use_t_plotting = [];
clearvars motion_energy_crop motion_energy motion_energy_std_crop motion_energy_std
[motion_energy, motion_energy_std, cam_align_diff_avg, n_trials] = cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

figure();
plot(surround_t(1:end-1), motion_energy); hold on;
plotshaded(surround_t(1:end-1), [motion_energy' - motion_energy_std'; ...
motion_energy' + motion_energy_std'], 'b');
ymax = ylim;
line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
xlabel('time from wheel movement onset')
ylabel('motion energy (whole image)')
legend(['# trials = ' num2str(n_trials)])
makepretty;

figure();
use_t_plotting = [0, 1];
use_t_idx = surround_t >= use_t_plotting(1) & surround_t <= use_t_plotting(2);   
imagesc(nanmean(cam_align_diff_avg(:, :, use_t_idx(2:end)), 3));

for iCrop = 1:size(ROIs,1)
    box_crop = ROIs(iCrop,:);
    plot_me = 0;
    use_t_plotting = [0, 0.2];
    [motion_energy_crop(iCrop,:), motion_energy_std_crop(iCrop,:), ~, n_trials_crop(iCrop)] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);
end

wheel_window = [-5, 5];
wheel_window_t = wheel_window(1):1 / framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus, use_align, wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, wheel_window_t_peri_event);
event_aligned_wheel_ds = resample(event_aligned_wheel, 10, 3);



axis image off;
%title(stim_types{iCondition});
paw_ROI = drawrectangle('Label','paw');
snout_ROI = drawrectangle('Label','snout');
whisking_ROI = drawrectangle('Label','whisking');
ROIs = [paw_ROI.Position; snout_ROI.Position; whisking_ROI.Position];

save([animal,'_' day, '_manualROIs'],"ROIs")

%% Get average movie around event
clear all;

animal = 'JF093';
day = '2023-03-06';
% go no go example: animal JF093, day 2023-03-06 exp 2
% go go go example: JF096, day 2023-02-07 exp 2
experiment = 2;
load_parts.cam = true;
load_parts.ephys = true;
verbose = true;
site = 1;
JF_load_experiment;
camera_location = facecam_fn;
camera_t = facecam_t;
% Get wheel movements during stim, only use quiescent trials

framerate = 30;

wheel_window = [0, 0.5];
wheel_window_t = wheel_window(1):1 / framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus, stimOn_times, wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0, 2);

timeline_rate = 0.0010;

stims = [4,12,6];
stim_types = {'Stim 1', 'Stim 2', 'Stim 3'};
%surround_frames1 = 30;
surround_frames = 150;
surround_t = [-surround_frames:surround_frames]./framerate;

%use_stim = stims(iStim);
use_align = stimOn_times(~isnan(stim_to_move)) + stim_to_move(~isnan(stim_to_move));%stimOn_times(stimIDs == use_stim & quiescent_trials & no_move_trials==1);

% Initialize video reader, get average and average difference
box_crop = [];
plot_me = 0;
use_t_plotting = [];
[motion_energy, motion_energy_std, cam_align_diff_avg, n_trials] = cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

figure();
plot(surround_t(1:end-1), motion_energy); hold on;
plotshaded(surround_t(1:end-1), [motion_energy' - motion_energy_std'; ...
motion_energy' + motion_energy_std'], 'b');
ymax = ylim;
line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
xlabel('time from wheel movement onset')
ylabel('motion energy (whole image)')
legend(['# trials = ' num2str(n_trials)])
makepretty;


%% sanity check basic 
% whole frame + "paw" region 

%use_stim = stims(iStim);
use_align = stimOn_times(~isnan(stim_to_move));% + stim_to_move(~isnan(stim_to_move));%stimOn_times(stimIDs == use_stim & quiescent_trials & no_move_trials==1);
use_align = use_align(400);
% Initialize video reader, get average and average difference
box_crop = [];
plot_me = 0;
use_t_plotting = [];
clearvars motion_energy_crop motion_energy motion_energy_std_crop motion_energy_std
[motion_energy, motion_energy_std, cam_align_diff_avg, n_trials] = cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);
for iCrop = 1:size(ROIs,1)
    box_crop = ROIs(iCrop,:);
    plot_me = 0;
    use_t_plotting = [0, 0.2];
    [motion_energy_crop(iCrop,:), motion_energy_std_crop(iCrop,:), ~, n_trials_crop(iCrop)] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);
end

wheel_window = [-5, 5];
wheel_window_t = wheel_window(1):1 / framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus, use_align, wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, wheel_window_t_peri_event);
event_aligned_wheel_ds = resample(event_aligned_wheel, 10, 3);

figure();
subplot(311); hold on;
s_idx = find(stimOn_times >= use_align - 5 & stimOn_times <= use_align + 5);
these_stims = stimOn_times(s_idx);
arrayfun(@(x) plot([these_stims(x) - use_align, these_stims(x) - use_align+0.4], [1,1]), 1:size(these_stims , 1));
ylabel('stimOn')

subplot(312)
%plot(surround_t(1:end-1), motion_energy); hold on;
%yyaxis left;
plot(surround_t(1:end-1), squeeze(motion_energy_crop(1,:))); 
ylabel('paw motion energy')
%yyaxis right;

subplot(313)
plot(-5:10/301:5-10/301, event_aligned_wheel)
ylabel('wheel movement')


%% sanity check: motion energy aligned to wheel movement onset
% whole frame + "paw" region 

%use_stim = stims(iStim);
use_align = stimOn_times(~isnan(stim_to_move)) + stim_to_move(~isnan(stim_to_move));%stimOn_times(stimIDs == use_stim & quiescent_trials & no_move_trials==1);

% Initialize video reader, get average and average difference
box_crop = [];
plot_me = 0;
use_t_plotting = [];
[motion_energy, motion_energy_std, cam_align_diff_avg, n_trials] = cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

figure();
plot(surround_t(1:end-1), motion_energy); hold on;
plotshaded(surround_t(1:end-1), [motion_energy' - motion_energy_std'; ...
motion_energy' + motion_energy_std'], 'b');
ymax = ylim;
line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
xlabel('time from wheel movement onset')
ylabel('motion energy (whole image)')
legend(['# trials = ' num2str(n_trials)])
makepretty;

%% define paw, snout, whisking regions 
figure();
 use_t_plotting = [0, 0.2];
        use_t_idx = surround_t >= use_t_plotting(1) & surround_t <= use_t_plotting(2);
      
imagesc(nanmean(cam_align_diff_avg(:, :, use_t_idx(2:end)), 3));

axis image off;
%title(stim_types{iCondition});
paw_ROI = drawrectangle('Label','paw');
snout_ROI = drawrectangle('Label','snout');
whisking_ROI = drawrectangle('Label','whisking');
ROIs = [paw_ROI.Position; snout_ROI.Position; whisking_ROI.Position];

save([animal,'_' day, '_manualROIs'],"ROIs")

for iCrop = 1:size(ROIs,1)
    box_crop = ROIs(iCrop,:);
    plot_me = 0;
    use_t_plotting = [0, 0.2];
    [motion_energy_crop(iCrop,:), motion_energy_std_crop(iCrop,:), ~, n_trials_crop(iCrop)] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);
end

%get ROIs mapped onto the frame pixels 
conds = {'whole video', 'paw', 'snout', 'whisking'};

figure();
subplot(1,4,1);
plot(surround_t(1:end-1), motion_energy); hold on;
plotshaded(surround_t(1:end-1), [motion_energy' - motion_energy_std'; ...
motion_energy' + motion_energy_std'], 'b');
ymax = ylim;
line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
legend(conds{1})
makepretty;
xlabel('time from wheel movement onset')
ylabel('motion energy')

colorsCrop = [rgb('DarkOrange'); rgb('HotPink'); rgb('DarkSlateGrey')];
for iCondition =1:3
    subplot(1,4,1+iCondition);
    plot(surround_t(1:end-1), motion_energy_crop(iCondition,:), 'Color', colorsCrop(iCondition,:));
    plotshaded(surround_t(1:end-1), [motion_energy_crop(iCondition,:) - motion_energy_std_crop(iCondition,:); ...
    motion_energy_crop(iCondition,:) + motion_energy_std_crop(iCondition,:)], colorsCrop(iCondition,:));
    ymax = ylim;
    line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
    legend(conds{1+iCondition})
    makepretty;
end

save([animal,'_' day, '_manualROIs_moveAlign_MI'],"motion_energy_crop", "motion_energy_std_crop")

%% aligned to stimulus onset 
for iStim = 1:3
    use_stim = stims(iStim);
    for iMove= 1:3
        if iMove == 1
            use_align = stimOn_times(stimIDs == use_stim & no_move_trials==1);
        elseif iMove==2
            use_align = stimOn_times(stimIDs == use_stim & no_move_trials==0);
        else
            use_align = stimOn_times(stimIDs == use_stim );
        end
    surround_frames = 90;
    surround_t = [-surround_frames:surround_frames]./framerate;
for iCrop = 1:size(ROIs,1)
    box_crop = ROIs(iCrop,:);
    plot_me = 0;
    use_t_plotting = [0, 0.2];
    [motion_energy_stim_crop(iStim, iCrop,iMove, :), motion_energy_std_stim_crop(iStim, iCrop, iMove, :), ~, n_trials_stim_crop(iStim, iCrop)] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);
end
    end
end

%get ROIs mapped onto the frame pixels 
conds = {'whole video', 'paw', 'snout', 'whisking'};

for iMove = 1:3
figure();
%iMove = 3;
for iStim = 1:3
    use_stim = stims(iStim);
% subplot(1,4,1);
% plot(surround_t(1:end-1), motion_energy); hold on;
% plotshaded(surround_t(1:end-1), [motion_energy' - motion_energy_std'; ...
% motion_energy' + motion_energy_std'], 'b');
% ymax = ylim;
% line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
% legend(conds{1})
% makepretty;
% xlabel('time from wheel movement onset')
% ylabel('motion energy')

colorsCrop = [rgb('Navy'); rgb('Crimson'); rgb('DarkSlateGrey')];

for iCondition =1:3
    
    subplot(1,3,iCondition); hold on;
    plot(surround_t(1:end-1), squeeze(motion_energy_stim_crop(iStim, iCondition, iMove, :)), 'Color', colorsCrop(iStim ,:));
    plotshaded(surround_t(1:end-1), [squeeze(motion_energy_stim_crop(iStim, iCondition, iMove, :))' - squeeze(motion_energy_std_stim_crop(iStim, iCondition, iMove, :))'; ...
    squeeze(motion_energy_stim_crop(iStim, iCondition, iMove, :))' + squeeze(motion_energy_std_stim_crop(iStim, iCondition, iMove, :))'], colorsCrop(iStim,:));
    ymax = ylim;
    line([0, 0],[ymax(1), ymax(2)], 'Color', 'red')
    title(conds{1+iCondition})
    %legend(conds{1+iCondition})
    makepretty;
end
end
end
save([animal,'_' day, '_manualROIs_MI'],"motion_energy_stim_crop", "motion_energy_std_stim_crop")

%% movement x visual response 

%% 1. get motion energy for every trial 
surround_frames = 30;
surround_t = [-surround_frames:surround_frames]./framerate;

% Initialize video reader, get average and average difference
box_crop = [];
plot_me = 0;
use_t_plotting = [];

stims = [4,12,6];
stimOn_times_use = stimOn_times(ismember(trial_conditions(:,1), stims)& ismember(trial_conditions(:,2), -90));
re_idx = ismember(trial_conditions(ismember(trial_conditions(:,1), stims),2), -90);
%rere = ismember(trial_conditions(:,1), stims);

    box_crop = [];
for iTrial = 1:size(stimOn_times_use,1)
    use_align = stimOn_times_use(iTrial);
     [motion_energy_stim_trial(iTrial, :), ~, ~, ~] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

end


for iCrop = 1:size(ROIs,1)
    box_crop = ROIs(iCrop,:);
for iTrial = 1:size(stimOn_times_use,1)
    use_align = stimOn_times_use(iTrial);
     [motion_energy_stim_trial_crop(iTrial, iCrop, :), ~, ~, ~] = ...
        cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting);

end
end
figure();
hold on;
for ii=1:10
    plot(squeeze(motion_energy_stim_trial(ii,:)))
end
%% 2. get average visual response for every trial 
all_temps = unique(spike_templates(spike_depths >=1600));
raster_window = [-1, 1];
clearvars ucrr_raster
psth_bin_size = 0.001;
for iTrial = 1:size(stimOn_times_use,1)
    [~,curr_raster(iTrial, :), t, ~, ~] = cl_raster_psth(spike_templates, spike_times_timeline, ...
        all_temps, raster_window, psth_bin_size, stimOn_times_use(iTrial), []);
end

%% 3. 
% bin movement [0, 0.3s] after stimulus onset. histogram 
figure(2); clf;
figure(3); clf;
crops = {'paw', 'snout', 'whisking'};
for iCrop = 1:3
    average_MI = nanmean(motion_energy_stim_trial_crop(:, iCrop, 30:40),3);
figure(2);
subplot(3,1,iCrop)
histogram(average_MI, 50)
xlabel('mean motion energy [0-0.3]s after stim onset')
ylabel('# trials')
title(crops{iCrop})
makepretty;

bins = [min(average_MI),  prctile(average_MI,20), prctile(average_MI,40), ...
    prctile(average_MI,60), prctile(average_MI,80),  max(average_MI)];
bins_center = nanmean([bins(1:end-1); bins(2:end)]);
[N,X] = hist(average_MI, bins);

Y = discretize(average_MI, bins);

%stimOn_times_use = stimOn_times(ismember(stimIDs, stims));
stimOn_ID_use = stimIDs(ismember(stimIDs, stims));
stims = [4,12,6];
for iBin = 1:max(Y)
    for iStim = 1:3
        %bin_average_mi(iStim, iBin) = nanmean(average_MI(Y==iBin & stimOn_ID_use==stims(iStim)));
        bin_average_vis(iStim, iBin) = nanmean(nanmean(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim) & re_idx, 1050:1150),2))*(1/0.15);
          bin_average_vis_std(iStim, iBin) = nanstd(nanmean(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim) & re_idx, 1050:1150),2)*(1/0.15)) ./...
            sqrt(size(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim) &re_idx,:),1));
    end
end

figure(3);
subplot(3,1,iCrop)
hold on;
for iStim = 1:3
    %subplot(1,3,iStim)
    plot(bins_center, squeeze(bin_average_vis(iStim,:)))
        plotshaded(bins_center, [squeeze(bin_average_vis(iStim,:)) - squeeze(bin_average_vis_std(iStim,:)); ...
        squeeze(bin_average_vis(iStim,:)) + squeeze(bin_average_vis_std(iStim,:))], cols(iStim,:))

end
    xlabel('motion energy')
    ylabel('visual response')
    makepretty;
legend({'go1', 'go2', 'no go'})
end

figure(4); clf;
figure(5); clf; 


average_MI = nanmean(motion_energy_stim_trial(:, 30:40),2);
figure(4);

histogram(average_MI, 50)
xlabel('mean motion energy [0-0.3]s after stim onset')
ylabel('# trials')

makepretty;

bins = [min(average_MI),  prctile(average_MI,20), prctile(average_MI,40), ...
    prctile(average_MI,60), prctile(average_MI,80),  max(average_MI)];
bins_center = nanmean([bins(1:end-1); bins(2:end)]);

[N,X] = hist(average_MI, bins);

Y = discretize(average_MI, bins);

%stimOn_times_use = stimOn_times(ismember(stimIDs, stims));
stimOn_ID_use = stimIDs(ismember(stimIDs, stims));
stims = [4,12,6];
clearvars bin_average_mi bin_average_vis

for iBin = 1:max(Y)
    for iStim = 1:3
        bin_average_mi(iStim, iBin) = nanmean(average_MI(Y==iBin & stimOn_ID_use==stims(iStim) &re_idx));
        bin_average_vis(iStim, iBin) = nanmean(nanmean(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim)&re_idx, 1050:1150),2))*(1/0.15);
        bin_average_vis_std(iStim, iBin) = nanstd(nanmean(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim)&re_idx, 1050:1150),2)*(1/0.15)) ./...
            sqrt(size(curr_raster(Y==iBin & stimOn_ID_use==stims(iStim)&re_idx,:),1));
    end
end

figure(5);clf;

hold on;
cols = [rgb('blue'); rgb('red'); rgb('Gold')];
for iStim = 1:3
    %subplot(1,3,iStim)
    sm = smoothdata(squeeze(bin_average_vis(iStim,:)));
    plot(bins_center, squeeze(bin_average_vis(iStim,:)))
    plotshaded(bins_center, [squeeze(bin_average_vis(iStim,:)) - squeeze(bin_average_vis_std(iStim,:)); ...
        squeeze(bin_average_vis(iStim,:)) + squeeze(bin_average_vis_std(iStim,:))], cols(iStim,:))
    

end
    xlabel('motion energy')
    ylabel('visual response')
    makepretty;
legend({'go1', '','go2', '','no go'})






