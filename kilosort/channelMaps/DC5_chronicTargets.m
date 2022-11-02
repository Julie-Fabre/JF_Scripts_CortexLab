
% plot just the best guys
[~, brain_outline] = plotBrainGrid([], []);
theseColors = {rgb('DeepSkyBlue');   rgb('DarkOrange'); rgb('Crimson')};
regionsNames = {'CP', 'GPe', 'SNr'};
for iRegion = 1:size(regionsNames,2)
    curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
 structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0); 
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
            'FaceColor', theseColors{iRegion}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
       
end
colsThese = lines(length(probe_ccf));
if size(colsThese,1) > 7
%    hashes = qq;
else
    
end
hold on;
plot3([480, 480],[690, 690],[48 621], 'Color', rgb('DeepSkyBlue'), 'LineWidth', 2) %str; 
plot3([630, 630],[810, 810],[48 621], 'Color', rgb('DarkOrange'), 'LineWidth', 2) % gpe
plot3([860, 860],[700, 700],[48 621], 'Color', rgb('Crimson'), 'LineWidth', 2) % snr


view([90, 0])

view([0, 0])

scatter3(480, 690, 48 , 120,  rgb('DeepSkyBlue'), 'filled', 'MarkerEdgeColor', 'k') %str; 
scatter3( 630, 810, 48 , 120,  rgb('DarkOrange'), 'filled','MarkerEdgeColor', 'k') % gpe
scatter3( 860,700, 48 ,120,   rgb('Crimson'), 'filled', 'MarkerEdgeColor', 'k') % snr
view([0, 90])