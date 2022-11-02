function [curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group, sort_by, plot_me)
% Example:
% ---------
% thisTemplate = 1;
% raster_window = [-0.5, 2];
% align_times = stimOn_times;
% align_group = stimIDs;
% psth_bin_size = 0.001;
% sort_by = stim_to_move;
% plot_me = true;
% [curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
%     thisTemplate, raster_window, psth_bin_size, align_times, align_group,...
    sort_by, plot_me)
% Set default raster times
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins) ./ 2;
use_align = reshape(align_times, [], 1);
t_peri_event = use_align + t_bins;
% (handle NaNs by setting rows with NaN times to 0)
t_peri_event(any(isnan(t_peri_event), 2), :) = 0;

curr_spikes_idx = ismember(spike_templates, thisTemplate);
curr_raster_spike_times = spike_times_timeline(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
    curr_raster_spike_times > max(t_peri_event(:))) = [];

if ~any(diff(reshape(t_peri_event', [], 1)) < 0)
    % (if no backward time jumps, can do long bin and cut out in-between, faster)
    curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
        reshape(t_peri_event', [], 1)), NaN], size(t_peri_event'))';
    curr_raster = curr_raster_continuous(:, 1:end-1);
else
    % (otherwise, bin trial-by-trial)
    curr_raster = cell2mat(arrayfun(@(x) ...
        histcounts(curr_raster_spike_times, t_peri_event(x, :)), ...
        [1:size(t_peri_event, 1)]', 'uni', false));
end
smooth_size = 51;


if ~isempty(align_group)


    if ~isempty(sort_by)
        [~, sortedAlignId] = sort(align_group(1:size(sort_by, 1)));
        % is sortbyThis is not empty, sort by those times
        [sortedSecVar, sortSecId] = sort(sort_by(sortedAlignId));

        [raster_y, raster_x] = find(curr_raster(sortSecId(sortedAlignId), :));
        if plot_me
            figure();
            scatter(t(raster_x), raster_y, 2, 'filled')
            hold on;
            plot(sortedSecVar(sortedAlignId), 1:max(raster_y))
            set(gca, 'YDir', 'reverse')
        end
    else
        [~, sortedAlignId] = sort(align_group);
        [raster_y, raster_x] = find(curr_raster(sortedAlignId, :));
        if plot_me
            figure();
            scatter(t(raster_x), raster_y, 2, 'filled')
            set(gca, 'YDir', 'reverse')
        end

    end
    gw = gausswin(smooth_size, 3)';
    %gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
    %- uncomment this ^ if you want to preserve onset times.

    smWin = gw ./ sum(gw);
    bin_t = mean(diff(t_bins));

    curr_psth = grpstats(curr_raster, align_group(1:size(curr_raster, 1)), @(x) mean(x, 1));
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;


else
    [raster_y, raster_x] = find(curr_raster);
    if plot_me
        figure();
        scatter(t(raster_x), raster_y, 2, 'filled')
        set(gca, 'YDir', 'reverse')
    end

    gw = gausswin(smooth_size, 3)';
    gw(1:round(smooth_size/2)) = 0;
    smWin = gw ./ sum(gw);
    bin_t = mean(diff(t_bins));

    curr_psth = nanmean(curr_raster);
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;
end
end