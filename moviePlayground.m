%% Get average movie around event
animal = 'JF096';
day = '2023-02-07';
experiment = 1;
site = 1;
load_parts.cam = true;
verbose = true;
JF_load_experiment;
use_cam = facecam_fn;
use_t = facecam_t;
% Get wheel movements during stim, only use quiescent trials
unexpect_rew = signals_events.responseTimes(trial_conditions(:,5)==1);

[allRew, allRew_sorted] = sort([ signals_events.responseTimes(n_trials), unexpected_reward_times]);
allRew_idx =  [trial_conditions(:,5); 2*ones(size(unexpected_reward_times,2),1)];
allRew_idx_sorted = allRew_idx(allRew_sorted);
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

%use_stim = 3;

alignNum = [0, 1, 2];
alignString = {'rewards', 'reward omissions', 'unexpected rewards'};
for iAlign = 1

use_align = allRew(allRew_idx_sorted == alignNum(iAlign));
%use_align = stimOn_times(stimIDs == use_stim & quiescent_trials);
surround_frames = 30;
% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);
cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);
cam_align_diff_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2);
frame_t_offset = nan(size(use_align));
for curr_align = 1:length(use_align)
    
    % Find closest camera frame to timepoint
    [frame_t_offset(curr_align),curr_frame] = ...
        min(abs(use_align(curr_align) - use_t));
    if curr_frame > 30
    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    curr_clip = double(squeeze(read(vr,curr_surround_frames)));
    curr_clip_diff = abs(diff(curr_clip,[],3));
    cam_align_avg = cam_align_avg + curr_clip./length(use_align);
    cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff./length(use_align);
    AP_print_progress_fraction(curr_align,length(use_align));
    end
end
    surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
    AP_imscroll(cam_align_avg,surround_t)
    axis image;
    surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
    AP_imscroll(cam_align_diff_avg,surround_t(2:end))
    axis image;
    % Plot difference within window
    use_t = [0,0.2];
    use_t_idx = surround_t >= use_t(1) & surround_t <= use_t(2);
    figure;
    imagesc(nanmean(cam_align_diff_avg(:,:,use_t_idx(2:end)),3));
    axis image off;
    title(alignString{iAlign});
    clearvars cam_align_diff_avg use_t_idx surround_t frame_t_offset cam_align_avg
end

