
%% load naive visual recordings
corona = 0;
mice = {'AP080', 'AP082', 'JF017', 'JF019'}; % JF020 doesn't have histology yet. JF021 and JF022 to come, and many more
locations = {'CP', 'SNr', 'GPe', 'GPi', 'STN'};
protocols = {'JF_natural_images', 'JF_locations', 'JF_GratingsPassive'};
ephysData = struct;

allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);
thisCount = 1;
for iMouse = 1:size(mice, 2)
    %find experiments
    animal = mice{iMouse};
    protocol = 'JF'; %contains this - refine protocol later
    experiments = AP_find_experimentsJF(animal, protocol, 1);
    experiments = experiments([experiments.ephys]);
    % check if kilosorted
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)

    for iDate = 1:size(experiments, 1)
        corrDay = find(arrayfun(@(x) probe2ephys(x).day == iDate, 1:numel(probe2ephys)));
        for iCorrDay = 1:size(corrDay, 2)
            sites(iCorrDay) = probe2ephys(corrDay(iCorrDay)).site;
        end
        for iSite = 1:size(sites, 2)
            corrSite = find(arrayfun(@(x) probe2ephys(x).site == iSite, 1:numel(probe2ephys)));
            for iLocation = 1:size(locations, 2)
                probe = intersect(corrDay, corrSite);
                if ~isempty(probe)
                    this_ccf = probe_ccf(probe); %find if contains location of interest
                    if ~any(structfun(@isempty, this_ccf))
                        theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
                        theseLocationsInterest = contains(theseLocations, locations{iLocation});
                        theseDepths = this_ccf.probe_depths(theseLocationsInterest);
                        if ~isempty(find(theseDepths))
                            %load experiment
                            curr_day = probe2ephys(probe).day;
                            day = experiments(curr_day).day;
                            experimentThese = experiments(curr_day).experiment; % experiment number
                            %experiment = experiment(probe2ephys(probe).site);
                            verbose = false; % display load progress and some info figures
                            load_parts.cam = false;
                            load_parts.imaging = false;
                            load_parts.ephys = true;
                            site = probe2ephys(probe).site;
                            for iExperiment = 1:size(experimentThese, 2)
                                experiment = experimentThese(iExperiment);
                                AP_load_experimentJF;
                                theseTemplates = template_depths >= min(theseDepths) & template_depths <= max(theseDepths); %correct depth units
                                %find closest depth and get closest location
                                A = repmat(this_ccf.probe_depths, [1, length(template_depths)]);
                                [minValue, closestIndex] = min(abs(A-template_depths'));
                                closestLocation = this_ccf.trajectory_coords(closestIndex, :);
                                %raster general and for each trial condition
                                theseSpikes = ismember(spike_templates, theseTemplates);
                                %save info
                                ephysData(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
                                ephysData(thisCount).spike_templates = spike_templates(theseSpikes);
                                ephysData(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
                                ephysData(thisCount).templatesOI = theseTemplates;
                                ephysData(thisCount).template_depths = template_depths(theseTemplates);
                                ephysData(thisCount).template_wf = waveforms(theseTemplates, waveform_peak(~isnan(new_spike_idx))- ...
                                    4:waveform_peak(~isnan(new_spike_idx))+4);
                                ephysData(thisCount).template_location = closestLocation;
                                ephysData(thisCount).protocol = expDef;
                                ephysData(thisCount).animal = animal;
                                ephysData(thisCount).date = day;
                                ephysData(thisCount).site = probe;
                                ephysData(thisCount).experiment = experiment;
                                ephysData(thisCount).location = locations(iLocation);
                                thisCount = thisCount + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end


%% plot responses - Striatum

%find == 'CP'
striatumData = find(arrayfun(@(x) contains(ephysData(x).location, 'CP'), 1:numel(ephysData)));

%find striatum limits - allen atlas
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure}, 2, [])') ./ 255;

%get limits and bins
binSize = 100;
structure_3d_lims = permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]);
[r, c, v] = ind2sub(size(structure_3d_lims), find(structure_3d_lims));
x_lim_um = [min(r), max(r)] .* 100;
y_lim_um = [min(c), max(c)] .* 100;
z_lim_um = [min(v), max(v)] .* 100;
dataBinsX = x_lim_um(1):binSize:x_lim_um(2);
dataBinsY = y_lim_um(1):binSize:y_lim_um(2);
dataBinsZ = z_lim_um(1):binSize:z_lim_um(2);

%find which bins striatal templates belong to 
uniqueRecordings = arrayfun(@(x) [ephysData(x).animal, ephysData(x).date, num2str(ephysData(x).site)], 1:numel(ephysData),'UniformOutput', false);%unique date + site
unique(uniqueRecordings, 'rows')
%allTemplates = 