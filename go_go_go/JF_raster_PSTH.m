function [curr_smoothed_psth, curr_psth, raster_x, raster_y, rasterColor, alignedVector] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group, sort_by, color_by, plot_me, causal_smoothing, figN)
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
%    sort_by, plot_me)
% Set default raster times
alignedVector = [];
rasterColor = [] ;
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


if ~isempty(align_group)
    if ~isempty(sort_by)
        [sortedSecVar, sortSecId] = sort(sort_by);  
        align_group = align_group(1:size(sort_by, 1));
        [sortedAlignVar, sortedAlignId] = sort(align_group(sortSecId));
        % is sortbyThis is not empty, sort by those times
        %[sortedSecVar, sortSecId] = sort(sort_by(sortedAlignId));   
        rasterColor = sortSecId(sortedAlignId);

        [raster_y, raster_x] = find(curr_raster(rasterColor, :));
        if plot_me
            figure(figN);
            clf;
            subplot(3, 1, [2:3])
            hold on;
            theseGroups = unique(sortedAlignVar);
            for iGroup = 1:length(theseGroups)
                plotThisGroup = find(ismember(sortedAlignVar, theseGroups(iGroup)));
                plotThisGroup_y = ismember(raster_y, plotThisGroup)
                scatter(t(raster_x(plotThisGroup_y)), raster_y(plotThisGroup_y), 2, 'filled')
            end
            
            alignedVector = sortedSecVar(sortedAlignId);
            plot(alignedVector, 1:length(sortedSecVar))
            set(gca, 'YDir', 'reverse')
            xlabel('time (s)')
            ylabel('trial (sorted)')
            makepretty;
            xlim([t(1), t(end)])
        end

    smooth_size = 51;
    gw = gausswin(smooth_size, 3)';
    if causal_smoothing
        gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
    end


    smWin = gw ./ sum(gw);
    bin_t = mean(diff(t_bins));

    curr_psth = grpstats(curr_raster, align_group(1:size(curr_raster, 1)), @(x) mean(x, 1));
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;
    if plot_me
        subplot(3, 1, 1)
        plot(t, curr_smoothed_psth)
        %xlabel('time (s)')
        ylabel('sp/s')
        makepretty;
    end
    
    else
        groups_unique = unique(align_group);
        [sortedGroup, sortedAlignId] = sort(align_group);
        [raster_y, raster_x] = find(curr_raster(sortedAlignId, :));
        for iGroup = 1:length(groups_unique)
                these_grp_rows = find(sortedGroup == groups_unique(iGroup));
        end
        rasterColor = these_grp_rows;

        if plot_me
            subplot(3, 1, [2:3])
            these_col_plot = lines(length(groups_unique));
            for iGroup = 1:length(groups_unique)
                these_grp_rows = find(sortedGroup == groups_unique(iGroup));
                plot_me = find(ismember(raster_y, these_grp_rows));
                scatter(t(raster_x(plot_me)), raster_y(plot_me), 2, these_col_plot(iGroup, :), 'filled');
                hold on;
            end
            set(gca, 'YDir', 'reverse')
            xlabel('time (s)')
            ylabel('trial # (sorted)')
            makepretty;
        end
    end
    smooth_size = 51;
    gw = gausswin(smooth_size, 3)';
    if causal_smoothing
        gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
    end


    smWin = gw ./ sum(gw);
    bin_t = mean(diff(t_bins));

    curr_psth = grpstats(curr_raster, align_group(1:size(curr_raster, 1)), @(x) mean(x, 1));
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;
    if plot_me
        subplot(3, 1, 1)
        plot(t, curr_smoothed_psth)
        %xlabel('time (s)')
        ylabel('sp/s')
        makepretty;
    end


else
    [raster_y, raster_x] = find(curr_raster);
    rasterColor = ones(length(align_times),1);
    smooth_size = 51;
    gw = gausswin(smooth_size, 3)';
    if causal_smoothing
        gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
    end


    smWin = gw ./ sum(gw);
    bin_t = mean(diff(t_bins));

    curr_psth = grpstats(curr_raster, ones(size(curr_raster, 1), 1), @(x) mean(x, 1));
    curr_smoothed_psth = conv2(padarray(curr_psth, ...
        [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
        smWin, 'valid') ./ bin_t;

    if plot_me
        clf;
        subplot(3, 1, [2:3])
        if ~isempty(color_by)
            theseCols = unique(color_by);
            these_col_plot = lines(length(theseCols));
            for iColor = 1:length(theseCols)
                these_col_rows = find(color_by == theseCols(iColor));
                plot_me = find(ismember(raster_y, these_col_rows));
                scatter(t(raster_x(plot_me)), raster_y(plot_me), 2, these_col_plot(iColor, :), 'filled');
                hold on;
            end
            set(gca, 'YDir', 'reverse')

            xlabel('time (s)')
            ylabel('unsorted trial #')
            makepretty;

            curr_psth = grpstats(curr_raster, color_by(1:size(curr_raster, 1)), @(x) mean(x, 1));
            curr_smoothed_psth = conv2(padarray(curr_psth, ...
                [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
                smWin, 'valid') ./ bin_t;

            subplot(3, 1, 1)
            plot(t, curr_smoothed_psth)
            xlabel('time (s)')
            ylabel('sp/s')
            makepretty;

        else
            scatter(t(raster_x), raster_y, 2, rgb('Green'),'filled')
            set(gca, 'YDir', 'reverse')
            gw = gausswin(smooth_size, 3)';
            if causal_smoothing
                gw(1:round(smooth_size/2)) = 0; %half gaussian to preserve onset times
            end
            xlabel('time (s)')
            ylabel('trial (unsorted)')
            makepretty;


            smWin = gw ./ sum(gw);
            bin_t = mean(diff(t_bins));

            curr_psth = nanmean(curr_raster);
            curr_smoothed_psth = conv2(padarray(curr_psth, ...
                [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
                smWin, 'valid') ./ bin_t;

            subplot(3, 1, 1)
            plot(t, curr_smoothed_psth, 'Color', rgb('Green'))
            xlabel('time (s)')
            ylabel('sp/s')
            makepretty;

        end


    end


end
end