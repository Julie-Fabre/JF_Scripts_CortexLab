
clear all;

animal = 'AP083';
day = '2020-10-19';
experiment = 1;
lilrig_load_experimentJF;
lilrig_retinotopyJF;
set(gcf,'color','w');
saveas(gcf, ['//znas.cortexlab.net/Subjects' filesep animal filesep day filesep 'retinotopy_' animal '_' day])
%AP_reference_outline('grid')

animal = 'AP084';
day = '2020-10-19';
experiment = 1;
lilrig_load_experimentJF;
lilrig_retinotopyJF;
set(gcf,'color','w');
saveas(gcf, ['//znas.cortexlab.net/Subjects' filesep animal filesep day filesep 'retinotopy_' animal '_' day])
%AP_reference_outline('grid')