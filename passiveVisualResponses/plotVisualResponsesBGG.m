keep ephysData

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
        plotBrainGrid();
        hold on;
        scatter3(bregma(1), bregma(3), bregma(2), 'filled')
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
            'FaceColor', theseColors{iRegion, :}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
      
       
        % 2-D
        ii = permute(av(1:slice_spacing:end, ...
            1:slice_spacing:end, 1:slice_spacing:end/2) == curr_plot_structure, [3, 1, 2]); % / 2 to only get one hemispehere
        [r, c, v] = ind2sub(size(ii), find(ii));
        
        figure();
        % x-y
        subplot(3,size(regions,2),iRegion)
        bb = boundary(r, c);
        plot(r(bb)-(bregma(1)/10), c(bb)-(bregma(3)/10), 'Color' ,theseColors{iRegion});
        xlabel('ML')
        ylabel('AP')
        axis equal
        axis square
        xyXLim = xlim;
        xyYLim = ylim;
        makepretty;

        % x-z
        subplot(3,size(regions,2),size(regions,2)+iRegion)
        bb = boundary(r, v);
        plot(r(bb)-(bregma(1)/10), -v(bb), 'Color' ,theseColors{iRegion});
        xlabel('ML')
        ylabel('DV')
        axis equal
        axis square
        makepretty;

        % y-z
        subplot(3,size(regions,2),size(regions,2)*2+(iRegion))
        bb = boundary(c, v);
        plot(c(bb)-(bregma(3)/10), -v(bb), 'Color' ,theseColors{iRegion});
        xlabel('AP')
        ylabel('DV')
        axis equal
        axis square
        makepretty;
        
        % bin region in 200um bins 
        xy = [r(bb), c(bb)];
        xyCountBins = {xyXLim(1)*10+bregma(3):100:xyXLim(2)*10+bregma(3), xyYLim(1)*10+bregma(1):100:xyYLim(2)*10+bregma(1)};
        
        % load all region data [-0.1 - 0.3]s after stim onset, with which bin belongs to 
        
        CNew = cat(1,ephysData.location);
        
        ind=find(ismember(CNew,'CP'));
        
        for iLocation = 1:size(ind,1)
            theseLocationsInfo(iLocation) = size(ephysData(ind(iLocation)).template_location,1);
        end
        theseLocationsInfo = cumsum(theseLocationsInfo);
        theseLocations = cat(1,ephysData(ind).template_location);
        theseLocationsBregmaAbs = [abs(theseLocations(:,3)-bregma(3)), abs(theseLocations(:,1)-bregma(1)), theseLocations(:,2)];
        %
        theseSpikeTimes = cat(1,ephysData(ind).spike_times_timeline);
        theseSpikeTemplates = cat(1,ephysData(ind).spike_templates);
        theseStimOnTimes = cat(ephysData(ind).stimOn_times,1);
        hold on;
        plot3(theseLocations(:,1),theseLocations(:,3),theseLocations(:,2),'.k')%1,2,3;
        [N,Xedges,Yedges,binX,binY] = histcounts2(theseLocations(:,1), theseLocations(:,2), xyCountBins{1,1},xyCountBins{1,2});%par rapport a bregma! 
        
        for iBinX = 1:size(Xedges,2)
            for iBinY = 1:size(Yedges,2)
                theseNeurons = binX==iBinX & binY==iBinY;
                if ~isempty(find(theseNeurons))
                    [h,hh] = histc(find(theseNeurons), theseLocationsInfo);
                    uniqueRecs = unique(hh);
                    for iUniqueRecs = uniqueRecs %get psth per rec 
                    end
                    
                    
                end 
                
            end
        end
        % average 
        
        % overlay response 

    end
end