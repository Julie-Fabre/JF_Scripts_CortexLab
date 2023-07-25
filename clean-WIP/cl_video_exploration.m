

%% Get average movie around event
animal = 'JF093';
day = '2023-03-06';
experiment = 2;
load_parts.cam = true;
verbose = true;
JF_load_experiment
use_cam = facecam_fn;
use_t = facecam_t;
% Get wheel movements during stim, only use quiescent trials
framerate = 30;
wheel_window = [0, 0.5];
wheel_window_t = wheel_window(1):1 / framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus, stimOn_times, wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity, wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0, 2);
stims = [4,12,6];
stim_types = {'Stim 1', 'Stim 2', 'Stim 3'};

%% general movements 

% snout ROI

figure();
plot( surround_t(use_t_store_idx(2:end)), motion_energy(:,1))
hold on;
plotshaded(surround_t(use_t_store_idx(2:end)), [motion_energy(:,1)' - motion_energy_std(:,1)'; ...
    motion_energy(:,1)' + motion_energy_std(:,1)'], 'b');
plot( surround_t(use_t_store_idx(2:end)), motion_energy(:,2))
plotshaded(surround_t(use_t_store_idx(2:end)), [motion_energy(:,2)' - motion_energy_std(:,2)'; ...
    motion_energy(:,2)' + motion_energy_std(:,2)'], 'r');
plot( surround_t(use_t_store_idx(2:end)), motion_energy(:,3)+nanmean(motion_energy(:,1))- ...
    nanmean(motion_energy(:,2)))
xx= motion_energy(:,3)' +nanmean(motion_energy(:,1))- ...
    nanmean(motion_energy(:,2)) - motion_energy_std(:,3)';
xxx = motion_energy(:,3)' +nanmean(motion_energy(:,1))- ...
    nanmean(motion_energy(:,2)) + motion_energy_std(:,3)';
plotshaded(surround_t(use_t_store_idx(2:end)), [xx; xxx], rgb('Gold'));
makepretty;

% whisker ROI

% loops 
for iStim = 1:3
    use_stim = stims(iStim);
    use_align = stimOn_times(stimIDs == use_stim & quiescent_trials);
    surround_frames = 30;
    % Initialize video reader, get average and average difference
    vr = VideoReader(use_cam);
    cam_im1 = read(vr, 1);
    cam_align_avg = zeros(size(cam_im1, 1), size(cam_im1, 2), ...
        surround_frames*2+1);
    cam_align_diff_avg = zeros(size(cam_im1, 1), size(cam_im1, 2), ...
        surround_frames*2);
    frame_t_offset = nan(size(use_align));
    for curr_align = 1:length(use_align)

        % Find closest camera frame to timepoint
        [frame_t_offset(curr_align), curr_frame] = ...
            min(abs(use_align(curr_align)-use_t));

        % Pull surrounding frames
        if curr_frame > 30
            curr_surround_frames = curr_frame + [-surround_frames, surround_frames];
            %curr_surround_frames = curr_surround_frames_seconds * framerate;
            curr_clip = double(squeeze(read(vr, curr_surround_frames)));
            curr_clip_diff = abs(diff(curr_clip, [], 3));
            cam_align_avg = cam_align_avg + curr_clip ./ length(use_align);
            cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff ./ length(use_align);
            use_t_storing = [-1, 1];
            use_t_store_idx = surround_t >= use_t_storing(1) & surround_t <= use_t_storing(2);
   
            motion_energy_frame(:,curr_align,iStim) = squeeze(nansum(nansum(abs(curr_clip_diff(:, :, use_t_store_idx(2:end))),1),2));
            AP_print_progress_fraction(curr_align, length(use_align));
        end
    end
    surround_t = [-surround_frames:surround_frames] ./ vr.FrameRate;
    AP_imscroll(cam_align_avg, surround_t)
    axis image;
    surround_t = [-surround_frames:surround_frames] ./ vr.FrameRate;
    AP_imscroll(cam_align_diff_avg, surround_t(2:end))
    axis image;
    % Plot difference within window
    use_t_plotting = [0, 0.2];
    use_t_idx = surround_t >= use_t_plotting(1) & surround_t <= use_t_plotting(2);
    figure;
    imagesc(nanmean(cam_align_diff_avg(:, :, use_t_idx(2:end)), 3));
    axis image off;
    title(stim_types{iStim});
    
    % store diff over time trace 
    use_t_storing = [-0.5, 0.5];
    use_t_store_idx = surround_t >= use_t_storing(1) & surround_t <= use_t_storing(2);
    motion_energy(:,iStim) = squeeze(nanmean(motion_energy_frame(:,:,iStim),2));
    motion_energy_std(:,iStim) = squeeze(nanstd(motion_energy_frame(:,:,iStim),[],2)./sqrt(size(motion_energy_frame,2)));
end


%motion_energy_frame(:,:,iStim))
%% snout movements 
figure;
imagesc(nanmean(cam_align_diff_avg(:, :, use_t_idx(2:end)), 3));
axis image off;
title(stim_types{iStim});
snout_ROI = drawrectangle;

%% whisking movements 

% roi = drawline