iProbe = iProbe+1
%for iProbe = 1:size(probe2ephys, 2)
keep st probe2ephys tv av animal iProbe slice_path im_path
    use_probe = iProbe;
    corona = 0;
    protocol = 'oice'; %protocol common to all sites and days
    experiments = AP_find_experimentsJF(animal, protocol, protocol);
    experiments = experiments([experiments.ephys]);
    curr_day = probe2ephys(iProbe).day;
    if isfield(probe2ephys, 'shank')
        curr_shank = probe2ephys(iProbe).shank;
    else
        curr_shank = NaN;
    end
    %if ~isempty(curr_day)
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment; % experiment number
    if length(experiment)>1
        experiment = experiment(2);
    end
    %experiment = experiment(probe2ephys(iProbe).site);
    verbose = false; % display load progress and some info figures
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    site = probe2ephys(iProbe).site;
    
    %experiment=[];
    
    lfp_channel = 'all'; % load lfp 
    loadClusters = 0;
    recording = [];
    [ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap',site,recording);
    isSpikeGlx = contains(ephysAPfile, '_g');%spike glx (2.0 probes) or open ephys (3A probes)? 
    if isSpikeGlx
         [ephysKSfile,~] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys',site,recording);
        if isempty(dir([ephysKSfile filesep 'sync.mat']))
            sync = syncFT(ephysAPfile, 385, ephysKSfile);
        end
        loadLFP=0;
    else
        loadLFP=1;
    end
    loadLFP=0;
    recording=[];
    experiment=1
    
    AP_load_experimentJF;

    %if dontAnalyze == 0
        AP_cellrasterJF({stimOn_times}, {stimIDs})
        
        
%         AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))',signals_events.responseTimes(n_trials(1):n_trials(end))'}, ...
%             {trial_conditions(:,1),trial_conditions(:,2), ...
%             trial_conditions(:,3),trial_outcome});
        lfp=nan;
        AP_align_probe_histologyJF(st, slice_path, ...
            spike_times, spike_templates, template_depths, spike_xdepths,template_xdepths,...
            lfp, lfp_channel_positions, lfp_channel_xpositions, ...
            use_probe,isSpikeGlx, curr_shank);