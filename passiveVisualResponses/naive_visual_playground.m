
%% plot striatum outline - "caudoputamen" from allen atlas
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure}, 2, [])') ./ 255;

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
structure_alpha = 0.2;

% 3-D plot
[~, brain_outline] = plotBrainGrid([], []);
hold on;
axis vis3d equal off manual
view([-30, 25]);
caxis([0, 300]);
[ap_max, dv_max, ml_max] = size(tv);
xlim([-10, ap_max + 10])
ylim([-10, ml_max + 10])
zlim([-10, dv_max + 10])
structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
scatter3(bregma(3), bregma(1), bregma(2));

% 2-D
ii = permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
[r, c, v] = ind2sub(size(ii), find(ii));

% x-y
figure();
bb = boundary(r, c);
plot(r(bb), c(bb));

% x-z
figure();
bb = boundary(r, v);
plot(r(bb), v(bb));

% y-z
figure();
bb = boundary(c, v);
plot(c(bb), v(bb));

%% visual projections
%get all V1 experiments
expIDs = findAllenExperiments('injection', 'VISp', 'line', '0', 'primary', true);
proj = getProjectionDataFromExperiment(expIDs(1));
%plot proj to caudoputamen, projection_intensity, projection volume,
%structure_id

[~, brain_outline] = plotBrainGrid([], []);
hold on;
axis vis3d equal off manual
view([-30, 25]);
caxis([0, 300]);
[ap_max, dv_max, ml_max] = size(tv);
xlim([-10, ap_max + 10])
ylim([-10, ml_max + 10])
zlim([-10, dv_max + 10])

%% GPi (AP083)
% load all data into one structure for gratings, location, nat images
animal = 'AP083';
AP_preprocess_phase3_newOEJF('AP083', '2020-12-05')
AP_preprocess_phase3_newOEJF('AP083', '2020-12-06')
AP_preprocess_phase3_newOEJF('AP083', '2020-12-07')

histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
load(histoFile)
probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
load(probe2ephysFile)
recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
recInfo = readtable(recInfoFile);

protocolStrings = {'rating', 'ocation', 'nat'}; 
%load probe
uniqueD = unique(recInfo.Date);
uniqueD = uniqueD(contains(uniqueD, '/')); %keep only dates (ie have a /) QQ
for iProbe = 1:size(probe2ephys, 2) %probe

    %get probe day and site. 
    curr_day = probe2ephys(iProbe).day; 
    curr_site = probe2ephys(iProbe).site; 
    day = uniqueD{curr_day};
    dayBreaks = find(day =='/'); 
    if dayBreaks(1)==2 && dayBreaks(2)==4
        day = strcat(day(dayBreaks(2)+1:end), '-0', day(dayBreaks(1)+1:dayBreaks(2)-1),...
        '-0', day(1:dayBreaks(1)-1)); %in standard format 
    elseif dayBreaks(1)==2
        day = strcat(day(dayBreaks(2)+1:end), '-', day(dayBreaks(1)+1:dayBreaks(2)-1),...
        '-0', day(1:dayBreaks(1)-1)); %in standard format 
    elseif dayBreaks(2)==5
        day = strcat(day(dayBreaks(2)+1:end), '-0', day(dayBreaks(1)+1:dayBreaks(2)-1),...
        '-', day(1:dayBreaks(1)-1)); %in standard format 
    else
        day = strcat(day(dayBreaks(2)+1:end), '-', day(dayBreaks(1)+1:dayBreaks(2)-1),...
        '-', day(1:dayBreaks(1)-1)); %in standard format 
    end
    
    theseExpDates = strcmp(recInfo.Date,uniqueD(curr_day));
    theseExpSites = recInfo.Site==curr_site;
    theseExpErrors = contains(recInfo.error, 'es') | ~isnan(recInfo.ephys_not_started) ...
        | ~isnan(recInfo.timeline_not_started); 
    theseExp = find(theseExpDates & theseExpSites & ~theseExpErrors); 
    theseProtocols = recInfo.Protocol(theseExp);
    %for each protocol type
    for iProtocol = 1:3
        experiment = recInfo.Protocol_number(theseExp(find(contains(theseProtocols, protocolStrings(iProtocol)))));
        
        AP_load_experimentJF;
        %
    end
end
