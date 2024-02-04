function [motion_energy, motion_energy_std, cam_align_diff_avg, n_trials,  motion_energy_frame_full] = cl_getMotionIndex(camera_location, camera_t, surround_frames, box_crop, surround_t, use_align, plot_me, use_t_plotting)
    % Initialize video reader, get average and average difference
    vr = VideoReader(camera_location);
    
   
    cam_im1 = read(vr, 1);
    % crop if necessary
    if ~isempty(box_crop)
        cam_im1 = cam_im1(round(box_crop(2)):round(box_crop(2)+box_crop(4)),...
                round(box_crop(1)):round(box_crop(1)+box_crop(3)));
    end
    cam_align_avg = zeros(size(cam_im1, 1), size(cam_im1, 2), ...
        surround_frames*2+1);
    cam_align_diff_avg = zeros(size(cam_im1, 1), size(cam_im1, 2), ...
        surround_frames*2);
    frame_t_offset = nan(size(use_align));

    for curr_align = 1:length(use_align)

        % Find closest camera frame to timepoint
        [frame_t_offset(curr_align), curr_frame] = ...
            min(abs(use_align(curr_align)-camera_t));

        % Pull surrounding frames
        if curr_frame > surround_frames
            curr_surround_frames = curr_frame + [-surround_frames, surround_frames];
            %curr_surround_frames = curr_surround_frames_seconds * framerate;
            curr_frames_clip = read(vr, curr_surround_frames);
            curr_clip = double(squeeze(curr_frames_clip));
            if ~isempty(box_crop)
                curr_clip = curr_clip(round(box_crop(2)):round(box_crop(2)+box_crop(4)),...
                    round(box_crop(1)):round(box_crop(1)+box_crop(3)),:);
            end
            curr_clip_diff = abs(diff(curr_clip, [], 3));
            cam_align_avg = cam_align_avg + curr_clip ./ length(use_align);
            cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff ./ length(use_align);
            use_t_storing = [surround_t(1), surround_t(end)];
            use_t_store_idx = surround_t >= use_t_storing(1) & surround_t <= use_t_storing(2);
   
            motion_energy_frame(:,curr_align) = squeeze(nansum(nansum(abs(curr_clip_diff(:, :, use_t_store_idx(2:end))),1),2));
            motion_energy_frame_full(:,curr_align) = squeeze(nansum(nansum(abs(curr_clip_diff(:, :, :)),1),2));
            %AP_print_progress_fraction(curr_align, length(use_align));
        end
    end
   
    if plot_me
        AP_imscroll(cam_align_avg, surround_t)
        title('Average aligned video frame')
        axis image;
        AP_imscroll(cam_align_diff_avg, surround_t(2:end))
        title('Average diff of aligned video frame')
        axis image;
        % Plot difference within window
        use_t_plotting = [0, 0.2];
        use_t_idx = surround_t >= use_t_plotting(1) & surround_t <= use_t_plotting(2);
        figure;
        imagesc(nanmean(cam_align_diff_avg(:, :, use_t_idx(2:end)), 3));
        axis image off;
    end
    
    % store diff over time trace 
    use_t_store_idx = surround_t >= use_t_storing(1) & surround_t <= use_t_storing(2);
    motion_energy = squeeze(nanmean(motion_energy_frame(:,:),2));
    motion_energy_std = squeeze(nanstd(motion_energy_frame(:,:),[],2)./sqrt(size(motion_energy_frame,2)));
    n_trials = size(motion_energy_frame,2);

end