animalsType = {'Naive'};
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
%regions = {'DMS', 'PS', 'STN', 'GPe', 'SNr', 'GPi'};
recordingInfo = readtable('C:\Users\Julie\Dropbox\Analysis\Recordings - Sheet1.csv');
bregma = [540,0,570];
%% X-Y, X-Z, Y-Z plots for each region. bin by 100 um (?) and plot in colormap increase in response ?

for iType = 1:size(animalsType, 2)
    theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    theseColors = {rgb('DeepSkyBlue'); rgb('SeaGreen'); rgb('DarkOrange'); rgb('Crimson'); rgb('Hotpink'); rgb('Black'); rgb('Brown')};

    allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
    tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
    av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
    slice_spacing = 10;
    structure_alpha = 0.2;
    figure();
    %overlay regions
    for iRegion = 1:size(regions, 2)

        curr_plot_structure = find(strcmp(st.acronym, regions{iRegion}));
        structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
        
        
       
        % 2-D
        ii = permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
        [r, c, v] = ind2sub(size(ii), find(ii));

        % x-y
        subplot(3,size(regions,2),iRegion)
        bb = boundary(r, c);
        plot(r(bb)-bregma(3)/10, c(bb)-bregma(1)/10, 'Color' ,theseColors{iRegion});
        xlabel('ML')
        ylabel('AP')
        axis equal
        axis square
        makepretty;

        % x-z
        subplot(3,size(regions,2),size(regions,2)+iRegion)
        bb = boundary(r, v);
        plot(r(bb)-bregma(3)/10, -v(bb), 'Color' ,theseColors{iRegion});
        xlabel('ML')
        ylabel('DV')
        axis equal
        axis square
        makepretty;

        % y-z
        subplot(3,size(regions,2),size(regions,2)*2+(iRegion))
        bb = boundary(c, v);
        plot(c(bb)-bregma(1)/10, -v(bb), 'Color' ,theseColors{iRegion});
        xlabel('AP')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
        
        % bin region in 200um bins 
        
        % load all region data [-0.1 - 0.3]s after stim onset, with which bin belongs to 

        % average 
        
        %overlay response 

    end
end