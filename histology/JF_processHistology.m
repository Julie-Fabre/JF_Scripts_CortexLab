%AP_preprocess_phase3_newOEJF_onlysync('JF022', '2020-12-14')

clear all;
animal = 'JF026';
animal = 'JF024';
im_type = 'brainSaw';

%% get histology slices and copy locally
if strcmp(im_type, 'brainSaw') == 1
     locationHisto = ['//znas.cortexlab.net/Brainsaw/JF026_JF029/', animal]; % copy files over to local disk
else
    locationHisto = ['//znas.cortexlab.net/Subjects/', animal, '/Histology/']; % copy files over to local disk
end
% Set paths for histology images and directory to save slice/alignment
if strcmp(im_type, 'o') == 1
    imageFolders = dir(locationHisto);
    imageFolders(1:2) = []; % remove current dir and one above
    im_dir = 'D:\';
    im_path = [im_dir, animal];
    slice_path = [im_dir, filesep, animal, filesep, 'slices'];
    mkdir(slice_path)
    for iFolder = 1:size(imageFolders, 1)
        copyfile([locationHisto, imageFolders(iFolder).name], im_path)
    end
elseif strcmp(im_type, 'brainSaw') == 1
    imagetFolders = dir([locationHisto, '/downsampled_stacks']);
    imagetFolders(1:2) = []; % remove current dir and one above
    imageFolders = dir([imagetFolders(1).folder '\' imagetFolders(1).name]);
    
    im_dir = 'D:\';
    im_path = [im_dir, animal];
    slice_path = [im_dir, filesep, animal, filesep, 'slices'];
%     mkdir(slice_path)
%     for iFolder = 1:size(imageFolders, 1)
%         imageFolders2 = dir([imageFolders(1).folder,filesep, '*.tif']);
%         for i = 1:size(imageFolders2,1)
%             copyfile([imageFolders(iFolder).folder, '\', imageFolders2(i).name], im_path)
%             theseImges = dir([im_path, '/img_*.tif']);
%             %movefile([im_path filesep theseImges.name],[im_path filesep 'img', num2str(1), 'channel', num2str(i), '.tif'],'f');
%         end
%     end
    
elseif strcmp(im_type, 'tiffUnmerged') == 1
    imageFolders = dir(locationHisto);
    imageFolders(1:2) = []; % remove current dir and one above

    %rename correctly 
    im_dir = 'D:\';
    im_path = [im_dir, animal];
    slice_path = [im_dir, filesep, animal, filesep, 'slices'];
    mkdir(slice_path)
    for iFolder = 1:size(imageFolders, 1)
        imageFolders2 = dir([locationHisto, imageFolders(iFolder).name, filesep, 'Default', filesep, '*.tif']);
        for i = 1:size(imageFolders2,1)
            copyfile([locationHisto, imageFolders(iFolder).name, filesep, 'Default', filesep, imageFolders2(i).name], im_path)
            theseImges = dir([im_path, '/img_*.tif']);
            movefile([im_path filesep theseImges.name],[im_path filesep 'img', num2str(iFolder), 'channel', num2str(i), '.tif'],'f');
        end
    end
elseif strcmp(im_type, 'tiffUnmergedNoDef') == 1
    imageFolders = dir(locationHisto);
    imageFolders(1:2) = []; % remove current dir and one above

    %rename correctly 
    im_dir = 'D:\';
    im_path = [im_dir, animal];
    slice_path = [im_dir, filesep, animal, filesep, 'slices'];
    mkdir(slice_path)
    for iFolder = 1:size(imageFolders, 1)
        imageFolders2 = dir([locationHisto, imageFolders(iFolder).name, filesep, '*.tif']);
        for i = 1:size(imageFolders2,1)
            copyfile([locationHisto, imageFolders(iFolder).name, filesep, imageFolders2(i).name], im_path)
            %theseImges = dir([im_path, '/img_*.tif']);
            %movefile([im_path filesep theseImges.name],[im_path filesep 'img', num2str(iFolder), 'channel', num2str(i), '.tif'],'f');
        end
    end
end

%% Allen atlas
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);

%% Set white balance and resize slide images, extract slice images
mkdir([im_path, filesep, 'resized'])
AP_process_histologyJF(im_path, im_type);

%% (optional) Rotate, pad, and center slice images
AP_rotate_histologyJF(slice_path);

%% align with ccf
AP_grab_histology_ccfJF(tv, av, st, slice_path);

% Align CCF slices and histology slices
% (first: automatically, by outline)
AP_auto_align_histology_ccfJF(slice_path,im_type);
% (second: curate manually)
AP_manual_align_histology_ccfJF(tv, av, st, slice_path);

%% 4) Utilize aligned CCF

% Display aligned CCF over histology slices
AP_view_aligned_histologyJF(st, slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volumeJF(tv, av, st, slice_path, 2);

% Get probe trajectory from histology, convert to CCF coordinates

AP_get_probe_histologyJF(tv, av, st, slice_path,'');

%% BLUE, ORANGE, YELLOW, PURPLE, GREEN , LIGHT BLUE, RED correspondance (manual):
probe2ephys = struct;
probe2ephys.animal = animal;
probe2ephys(1).day = 3;
probe2ephys(1).site = 2;

probe2ephys(2).day = 4;
probe2ephys(2).site = 2;

probe2ephys(3).day = 2;
probe2ephys(3).site = 1;

probe2ephys(4).day = 1;
probe2ephys(4).site = 1;

probe2ephys(5).day = 4;
probe2ephys(5).site = 1;

probe2ephys(6).day = 3;
probe2ephys(6).site = 1;

probe2ephys(7).day = 5;
probe2ephys(7).site = 1;

probe2ephys(8).day = 5;
probe2ephys(8).site = 2;

 probe2ephys(9).day = 3;
 probe2ephys(9).site = 2;
 
 probe2ephys(10).day = 4;
 probe2ephys(10).site = 2;

save([im_path, '/probe2ephys.mat'], 'probe2ephys')

%% align ephys
% Align histology to electrophysiology %% NOTE: YOU HAVE TO DO THESE ONE BY
% ONE RN (NOT IN LOOP) FOR THINGS TO SAVE PROPERLY
load([im_path, '/probe2ephys.mat'])
dontAnalyze = 0;
for iProbe = 1:size(probe2ephys, 2)
    use_probe = iProbe;
    corona = 0;
    protocol = 'vanilla'; %protocol common to all sites and days
    experiments = AP_find_experimentsJF(animal, protocol, protocol);
    experiments = experiments([experiments.ephys]);
    curr_day = probe2ephys(iProbe).day;
    if ~isempty(curr_day)
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
    
    lfp_channel = 'all';
    loadClusters = 0;

    isSpikeGlx=0;
    recording = 2;

    AP_load_experimentJF;

    if dontAnalyze == 0
        AP_cellrasterJF({stimOn_times}, {stimIDs})
        AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))'}, ...
    {trial_conditions(:,1).*trial_conditions(:,2), ...
    trial_conditions(:,3),trial_outcome});
        AP_align_probe_histologyJF(st, slice_path, ...
            spike_times, spike_templates, template_depths, ...
            lfp, channel_positions(:, 2), ...
            use_probe);
    else
        experiment = experiments(curr_day).experiment; % experiment number
        experiment = experiment(probe2ephys(iProbe).site) + 1;
        verbose = false; % display load progress and some info figures
        load_parts.cam = false;
        load_parts.imaging = false;
        load_parts.ephys = true;
        site = probe2ephys(iProbe).site;
        lfp_channel = 'all';
        try
            AP_load_experimentJF;
            AP_align_probe_histologyJF(st, slice_path, ...
                spike_times, spike_templates, template_depths, ...
                lfp, channel_positions(:, 2), ...
                use_probe);
        catch
            experiment = experiments(curr_day).experiment; % experiment number
            experiment = experiment(probe2ephys(iProbe).site) + 2;
            verbose = false; % display load progress and some info figures
            load_parts.cam = false;
            load_parts.imaging = false;
            load_parts.ephys = true;
            site = probe2ephys(iProbe).site;
            lfp_channel = 'all';
            AP_load_experimentJF;
            AP_align_probe_histologyJF(st, slice_path, ...
                spike_times, spike_templates, template_depths, ...
                lfp, channel_positions(:, 2), ...
                use_probe);
        end
    end
    end

    %     indexes = probe_ccf .* 10;
    %     av(indexes)
    %     av(round(probe_ccf(:, 1)), probe_ccf(:, 2), probe_ccf(:, 3))
end

%% save on server
mkdir([locationHisto, '/processed/']);
copyfile(im_path, [locationHisto, '/processed/']);