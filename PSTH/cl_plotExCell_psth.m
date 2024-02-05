function cl_plotExCell_psth(animal, experiment, iProbe, unit, align_type, group_type, ...
    raster_window, psth_bin_size, plot_raster, plot_psth, color_type)
% QQ to do: add s.e.m. , psths in spikes/per second 

cl_myPaths;

protocol = ''; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
experiments = experiments([experiments.ephys]);

orientationType = 'psl'; % psl (for posterior, superior, left), means the first, top left voxel
% is the most posterior, superior, left part of the brain
channelColToRegister = 'green'; % channel you want to register
channelColToTransform = 'red'; % channel you want to use to draw probes
atlasResolution_um = 25; % voxel size. currently in all dimensions, both for atlas and acquired image
atlasSpecies = 'mouse'; % atlas species
atlasType = 'allen'; % atlas name
brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives

brainsawPath_curr = [cl_cortexlab_filename(animal, '', '', 'histo_folder', '', '', ''), '/downsampled_stacks/025_micron'];
% registration location/files
[~, ~, ~, outputDir] = ...
    ya_getLocations(brainglobeLocation, brainsawPath_curr, channelColToRegister, ...
    channelColToTransform, atlasType, atlasSpecies, atlasResolution_um);
load([outputDir, '/probe2ephys.mat'])
load([outputDir, '/probe_ccf.mat'])

% load experiment
curr_day = probe2ephys(iProbe).day; % (set which day to use)
site = probe2ephys(iProbe).site;
curr_shank = probe2ephys(iProbe).shank;
if isfield(probe2ephys, 'recording')
    recording = probe2ephys(iProbe).recording;
else
    recording = [];
end
thisDate = experiments(curr_day).thisDate; % date
cl_load_experiment;

% get align times and group
if contains(align_type, 'stimOn_noMove')
    align_times = stimOn_times(no_move_trials);
    align_group = trial_conditions(no_move_trials, group_type);
elseif contains(align_type, 'stimOn')
    align_times = stimOn_times;
    align_group = trial_conditions(:, group_type);
end

% rename and adjust for shank # 
theseChannelPositions = [(curr_shank-1) * 250, (curr_shank-1)*250 + 32];
theseChannels = ismember(channel_positions(:,1), theseChannelPositions);
theseTemplates = ismember(template_xdepths, theseChannelPositions);
theseSpikes = ismember(spike_xdepths, theseChannelPositions);
spike_times_timeline = spike_times_timeline(theseSpikes);
spike_templates = spike_templates(theseSpikes);

good_templates_idx = unique(spike_templates);
new_spike_idx = nan(max(spike_templates), 1);
new_spike_idx(good_templates_idx) = 1:length(good_templates_idx);
spike_templates = new_spike_idx(spike_templates);

% get raster + psth
[curr_psth, ~, t, raster_x, raster_y, raster_group_id, psth_group_id] = cl_raster_psth(spike_templates, spike_times_timeline, ...
    unit, raster_window, psth_bin_size, align_times, align_group);

% plot
figure();
if plot_raster && plot_psth
    theseColors = cl_paper_plotting(color_type);
    subplot(5, 1, 1:4)
    raster_dots = scatter(t(raster_x), raster_y, 0.5, 'k');
    [~, mappedVector] = ismember(raster_group_id, psth_group_id);
    set(raster_dots,'CData',theseColors(mappedVector,:));
    xlim([raster_window(1)+psth_bin_size*10, raster_window(2)-psth_bin_size*10])

    subplot(5, 1, 5)
    for iGroup = 1:size(psth_group_id,1)
        plot(t,smoothdata(curr_psth(iGroup, :), 'gaussian', [0 50]), 'Color', theseColors(iGroup,:)); hold on;
    end
    xlim([raster_window(1)+psth_bin_size*10, raster_window(2)-psth_bin_size*10])
    hold off;

elseif plot_psth
    theseColors = cl_paper_plotting(color_type);
    for iGroup = 1:size(psth_group_id,1)
        plot(t,smoothdata(curr_psth(iGroup, :), 'gaussian', [0 50]), 'Color', theseColors(iGroup,:)); hold on;
    end
    xlim([raster_window(1)+psth_bin_size*10, raster_window(2)-psth_bin_size*10])
    hold off;
elseif plot_raster
    theseColors = cl_paper_plotting(color_type);
    subplot(5, 1, 1:4)
    raster_dots = scatter(t(raster_x), raster_y, 0.5, 'k');
    [~, mappedVector] = ismember(raster_group_id, psth_group_id);
    set(raster_dots,'CData',theseColors(mappedVector,:));
    xlim([raster_window(1)+psth_bin_size*10, raster_window(2)-psth_bin_size*10])

end


end