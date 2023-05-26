function [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group)


time_bins = raster_window(1):psth_bin_size:raster_window(2);
t = time_bins(1:end-1) + diff(time_bins) ./ 2;
use_align = reshape(align_times, [], 1);
t_peri_event = use_align + time_bins;
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
        (1:size(t_peri_event, 1))', 'uni', false));
end
if ~isempty(align_group)
    [~, sortedAlignId] = sort(align_group); 

    [raster_y, raster_x] = find(curr_raster(sortedAlignId, :));

    curr_psth = grpstats(curr_raster, align_group(1:size(curr_raster, 1)), @(x) mean(x, 1));
else
    curr_psth = mean(curr_raster);
    [raster_y, raster_x] = find(curr_raster(:, :));

end



end