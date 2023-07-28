function [motion_energy, motion_energy_std, n_trials] = cl_getMotionIndex(camera_location, camera_t, box_crop, surround_t, use_align)
 surround_frames = 30;
    % Initialize video reader, get average and average difference
    vr = VideoReader(camera_location);
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
   
            motion_energy_frame(:,curr_align,iCondition) = squeeze(nansum(nansum(abs(curr_clip_diff(:, :, use_t_store_idx(2:end))),1),2));
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
    title(stim_types{iCondition});
    
    % store diff over time trace 
    use_t_storing = [-0.5, 0.5];
    use_t_store_idx = surround_t >= use_t_storing(1) & surround_t <= use_t_storing(2);
    motion_energy(:,iCondition) = squeeze(nanmean(motion_energy_frame(:,:,iCondition),2));
    motion_energy_std(:,iCondition) = squeeze(nanstd(motion_energy_frame(:,:,iCondition),[],2)./sqrt(size(motion_energy_frame,2)));
    n_trials(iCondition) = size(motion_energy_frame,2);

end