

keep passive_data
%% log image
figure();
for iRegion = 1:3
     these_units = passive_data.unit_area==iRegion  &...
        (passive_data.unitType' ==1 | passive_data.unitType' ==2);

     subplot(3,3,(iRegion-1)*3 + 1)
     hist3([passive_data.templateDuration(these_units)', passive_data.pss(these_units)'],...
         'CDataMode','auto','FaceColor','interp','Edges', {150:50:1200 0:10:200})
     view(2);
     set(gca,'ColorScale','log')
     colorbar;
     xlabel('peak-to-trough duration (us)')
     ylabel('post spike suppression (ms)')

     subplot(3,3,(iRegion-1)*3 + 2)
     hist3([passive_data.templateDuration(these_units)', passive_data.propLongISI(these_units)'],...
         'CDataMode','auto','FaceColor','interp','Edges', {150:50:1200 0:0.05:1})
     view(2)
     colorbar;
     set(gca,'ColorScale','log')
     xlabel('peak-to-trough duration (us)')
     ylabel('frac. ISI > 2s')

     subplot(3,3,(iRegion-1)*3 + 3)
     hist3([passive_data.pss(these_units)', passive_data.propLongISI(these_units)'],...
         'CDataMode','auto','FaceColor','interp','Edges', {0:10:200 0:0.05:1})
     view(2)
     colorbar;
     set(gca,'ColorScale','log')
     xlabel('post spike suppression (ms)')
     ylabel('frac. ISI > 2s')

end
prettify_plot;

passive_data = load('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat');
passive_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
passive_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
%% a la Peters et al., 2021 
figure();
cellT_colors = ya_getColors(4);

for iRegion = 1:3
     these_units = passive_data.unit_area==iRegion  &...
        (passive_data.unitType' ==1 | passive_data.unitType' ==2);

     [regionClassification, unitClassification] = cl_subsection_region(passive_data.unit_coords,...
         passive_data.unit_area, passive_data.pss, passive_data.templateDuration, passive_data.propLongISI);
     unique_cellTypes = unique(unitClassification(passive_data.unit_area ==iRegion));
     unique_cellTypes(unique_cellTypes == "") =[];

     for iCellType = 1:size(unique_cellTypes,1)
     subplot(3,1,(iRegion-1)*1 + 1); hold on;
     scatter3(passive_data.templateDuration(these_units & unitClassification == unique_cellTypes(iCellType))',...
         passive_data.pss(these_units & unitClassification == unique_cellTypes(iCellType))', ...
         passive_data.propLongISI(these_units & unitClassification == unique_cellTypes(iCellType))', ...
         4000/sum(these_units &  unitClassification == unique_cellTypes(iCellType)),...
         cellT_colors(iCellType,:),'filled')
     xlim([150, 800])
     ylim([0 200])
     zlim([0, 1])

     xlabel('peak-to-trough duration (us)')
     ylabel('post spike suppression (ms)')
     zlabel('frac. ISI > 2s')
     set(gca, 'XDir', 'reverse')
     end

end
prettify_plot;
%xlim([150, 800])
%ylim([0 200])
%zlim([0, 1])
%% scatter
figure();
for iRegion = 1:3
     these_units = passive_data.unit_area==iRegion  &...
        (passive_data.unitType' ==1 | passive_data.unitType' ==2);

     subplot(3,3,(iRegion-1)*3 + 1)
     scatter(passive_data.templateDuration(these_units)', passive_data.pss(these_units)', 2, 'filled')


     subplot(3,3,(iRegion-1)*3 + 2)
     scatter(passive_data.templateDuration(these_units)', passive_data.propLongISI(these_units)', 2, 'filled')


     subplot(3,3,(iRegion-1)*3 + 3)
     scatter(passive_data.pss(these_units)', passive_data.propLongISI(these_units)', 2, 'filled')




end