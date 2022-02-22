
function [raster_x, raster_y, t, curr_smoothed_psth,trial_sort,curr_raster_sorted]= getRaster(spike_templates,thisUnit,spike_timeline, thisAlign, thisSort)
raster_window = [-0.3,0.5];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
t_peri_event = thisAlign + t_bins;
% (handle NaNs by setting rows with NaN times to 0)
t_peri_event(any(isnan(t_peri_event),2),:) = 0;

% Bin spikes (use only spikes within time range, big speed-up)
curr_spikes_idx = ismember(spike_templates,thisUnit);
curr_raster_spike_times = spike_timeline(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
    curr_raster_spike_times > max(t_peri_event(:))) = [];

if ~any(diff(reshape(t_peri_event',[],1)) < 0)
    % (if no backward time jumps, can do long bin and cut out in-between, faster)
    curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
        reshape(t_peri_event',[],1)),NaN],size(t_peri_event'))';
    curr_raster = curr_raster_continuous(:,1:end-1);   
else
    % (otherwise, bin trial-by-trial)
    curr_raster = cell2mat(arrayfun(@(x) ...
        histcounts(curr_raster_spike_times,t_peri_event(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false));
end

smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_t = mean(diff(t));

curr_psth =  mean(curr_raster,1);
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;
[~,trial_sort] = sort(thisSort);
curr_raster_sorted = curr_raster(trial_sort,:);
[raster_y,raster_x] = find(curr_raster_sorted);
end