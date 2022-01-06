function [curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, t_peri_event, timeT, stimGroup)
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
[raster_y, raster_x] = find(curr_raster);
if exist('stimGroup', 'var')
    gw = gausswin(smooth_size, 3)';
    %gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
    smWin = gw ./ sum(gw);
    bin_t = mean(diff(timeT));

    curr_psth = grpstats(curr_raster, stimGroup(1:size(curr_raster, 1)), @(x) mean(x, 1));
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;
  

else


    [~, sortMove] = sort(stim_to_move);
    [raster_y, raster_x] = find(curr_raster(sortMove, :));
    gw = gausswin(smooth_size, 3)';
    gw(1:round(smooth_size/2)) = 0;
    smWin = gw ./ sum(gw);
    bin_t = mean(diff(timeT));

    curr_psth = nanmean(curr_raster);
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;
end
end