function cl_plotExCell_psth(animal, experiment, iProbe, unit, align_type, group_type,...
    raster_window, psth_bin_size, plot_raster, plot_psth)

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
    recording =[];
end
thisDate = experiments(curr_day).thisDate; % date
cl_load_experiment;

% get align times and group 
if contains(align_type, 'stimOn_noMove')
    align_times = stimOn_times(no_move_trials);
    align_group = trial_conditions(no_move_trials,group_type);
elseif contains(align_type, 'stimOn')
    align_times = stimOn_times;
    align_group = trial_conditions(:,group_type);
end

% get raster + psth 
[curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
    iUnit, raster_window, psth_bin_size, align_times, align_group);

% plot 



end