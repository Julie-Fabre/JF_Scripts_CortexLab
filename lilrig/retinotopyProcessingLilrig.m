
clear all;

animal = 'JF023';
day = '2021-01-21';
experiment = 2;
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